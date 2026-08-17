#!/usr/bin/env python3
"""Generate the embedded C++ representation of the deployment XGBoost models."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import joblib


FEATURE_NAMES = [
    "mean_ref_ratio",
    "mean_alt_ratio",
    "mean_ref_count",
    "mean_alt_count",
    "ref_mix_alt_pure_frac",
    "ref_pure_alt_mix_frac",
    "both_mix_frac",
    "mean_abs_diff",
    "mean_abs_count_diff",
    "mean_center_dist_to_diag",
]

EXPECTED_THRESHOLDS = {
    "Snv": 0.79,
    "Indel": 0.17,
}

EXPECTED_MODEL_SHA256 = {
    "Snv": "65d8a81373640cffc4f29473a31407701ed8bd4816aea6d85cad5baf804de556",
    "Indel": "3adfd11447a22227b93dd7ae0bc229f585cf4a007e1d7eb536c1ccf2cb0b150c",
}

EXPECTED_METADATA_SHA256 = {
    "Snv": "bbb1ac41f3e894f5b28285ce59b10756d16dc7546ca9a0bc0a2d62efb52aba97",
    "Indel": "84e41b707297f938e4dd68b1174c7139eb86c534974747203b9a585bf8f2f3e6",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snv-model", required=True, type=Path)
    parser.add_argument("--indel-model", required=True, type=Path)
    parser.add_argument("--snv-metadata", type=Path)
    parser.add_argument("--indel-metadata", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str, description: str) -> str:
    actual = sha256(path)
    if actual != expected:
        raise ValueError(
            f"Unexpected {description} SHA-256 for {path}: "
            f"expected={expected}, actual={actual}"
        )
    return actual


def load_booster_json(path: Path) -> dict[str, Any]:
    model = joblib.load(path)
    if not hasattr(model, "get_booster"):
        raise TypeError(f"{path} is not an XGBoost sklearn estimator")

    actual_features = list(getattr(model, "feature_names_in_", []))
    if actual_features != FEATURE_NAMES:
        raise ValueError(
            f"Unexpected feature order in {path}:\n"
            f"expected={FEATURE_NAMES}\nactual={actual_features}"
        )

    return json.loads(model.get_booster().save_raw(raw_format="json"))


def load_metadata(path: Path, expected_threshold: float) -> dict[str, Any]:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if metadata.get("feature_names") != FEATURE_NAMES:
        raise ValueError(f"Unexpected metadata feature order in {path}")

    threshold = float(metadata.get("final_selected_threshold", math.nan))
    if threshold != expected_threshold:
        raise ValueError(
            f"Unexpected selected threshold in {path}: "
            f"expected={expected_threshold}, actual={threshold}"
        )
    return metadata


def parse_base_score(value: str) -> float:
    probability = float(value.strip("[]"))
    if not 0.0 < probability < 1.0:
        raise ValueError(f"Expected logistic base_score, got {value}")
    return math.log(probability / (1.0 - probability))


def float_literal(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError(f"Cannot embed non-finite float: {value}")
    literal = format(value, ".9g")
    if "e" not in literal.lower() and "." not in literal:
        literal += ".0"
    return literal + "f"


def double_literal(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError(f"Cannot embed non-finite double: {value}")
    literal = format(value, ".17g")
    if "e" not in literal.lower() and "." not in literal:
        literal += ".0"
    return literal


def render_model(symbol: str, model_json: dict[str, Any]) -> list[str]:
    learner = model_json["learner"]
    if learner["feature_names"] != FEATURE_NAMES:
        raise ValueError(f"Unexpected booster feature order for {symbol}")
    if learner["objective"]["name"] != "binary:logistic":
        raise ValueError(f"Unexpected objective for {symbol}")

    booster = learner["gradient_booster"]
    if booster["name"] != "gbtree":
        raise ValueError(f"Only gbtree boosters are supported, got {booster['name']}")

    model = booster["model"]
    trees = model["trees"]
    if any(tree_info != 0 for tree_info in model["tree_info"]):
        raise ValueError(f"Only binary single-output trees are supported for {symbol}")

    roots: list[int] = []
    rendered_nodes: list[str] = []
    node_offset = 0
    for tree in trees:
        roots.append(node_offset)
        left_children = tree["left_children"]
        right_children = tree["right_children"]
        split_indices = tree["split_indices"]
        split_conditions = tree["split_conditions"]
        default_left = tree["default_left"]
        split_type = tree["split_type"]

        if any(split_type):
            raise ValueError(f"Categorical splits are not supported for {symbol}")

        for node_idx, left_child in enumerate(left_children):
            right_child = right_children[node_idx]
            if left_child == -1:
                rendered_nodes.append(
                    "    {-1, 0.0f, -1, -1, -1, "
                    f"{float_literal(split_conditions[node_idx])}}},"
                )
                continue

            left_global = node_offset + left_child
            right_global = node_offset + right_child
            missing_global = left_global if default_left[node_idx] else right_global
            rendered_nodes.append(
                "    {"
                f"{split_indices[node_idx]}, "
                f"{float_literal(split_conditions[node_idx])}, "
                f"{left_global}, {right_global}, {missing_global}, 0.0f"
                "},"
            )

        node_offset += len(left_children)

    base_margin = parse_base_score(learner["learner_model_param"]["base_score"])
    lines = [
        f"const int k{symbol}Roots[{len(roots)}] = {{",
        "    " + ", ".join(str(root) for root in roots),
        "};",
        f"const MethylXgbNode k{symbol}Nodes[{len(rendered_nodes)}] = {{",
        *rendered_nodes,
        "};",
        f"const MethylXgbModel k{symbol}Model = {{",
        f"    k{symbol}Nodes, k{symbol}Roots, {len(roots)}, {double_literal(base_margin)}",
        "};",
    ]
    return lines


def main() -> None:
    args = parse_args()
    snv_metadata = args.snv_metadata or args.snv_model.with_name("final_selected_model_features.json")
    indel_metadata = args.indel_metadata or args.indel_model.with_name("final_selected_model_features.json")

    verify_sha256(args.snv_model, EXPECTED_MODEL_SHA256["Snv"], "Snv joblib")
    verify_sha256(args.indel_model, EXPECTED_MODEL_SHA256["Indel"], "Indel joblib")
    verify_sha256(snv_metadata, EXPECTED_METADATA_SHA256["Snv"], "Snv metadata")
    verify_sha256(indel_metadata, EXPECTED_METADATA_SHA256["Indel"], "Indel metadata")

    models = [
        ("Snv", args.snv_model, snv_metadata, load_booster_json(args.snv_model)),
        ("Indel", args.indel_model, indel_metadata, load_booster_json(args.indel_model)),
    ]

    for symbol, _, metadata_path, _ in models:
        load_metadata(metadata_path, EXPECTED_THRESHOLDS[symbol])

    lines = [
        "// Auto-generated by tools/generate_methyl_xgb_model_data.py.",
        "// Do not edit by hand.",
    ]
    for symbol, path, metadata_path, _ in models:
        lines.append(f"// {symbol} joblib SHA-256: {sha256(path)}")
        lines.append(f"// {symbol} metadata SHA-256: {sha256(metadata_path)}")
        lines.append(f"// {symbol} selected threshold: {EXPECTED_THRESHOLDS[symbol]}")
    lines.append("")

    for symbol, _, _, model_json in models:
        lines.extend(render_model(symbol, model_json))
        lines.append("")

    args.output.write_text("\n".join(lines), encoding="ascii")


if __name__ == "__main__":
    main()
