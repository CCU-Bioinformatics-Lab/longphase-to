#ifndef METHYL_XGB_MODEL_H
#define METHYL_XGB_MODEL_H

enum MethylXgbVariantKind {
    METHYL_XGB_KIND_SNV,
    METHYL_XGB_KIND_INDEL
};

// final_selected_model_features.json SHA-256:
// SNV   bbb1ac41f3e894f5b28285ce59b10756d16dc7546ca9a0bc0a2d62efb52aba97
// indel 84e41b707297f938e4dd68b1174c7139eb86c534974747203b9a585bf8f2f3e6
constexpr double METHYL_XGB_DEFAULT_SNV_THRESHOLD = 0.44;
constexpr double METHYL_XGB_DEFAULT_INDEL_THRESHOLD = 0.17;

struct MethylXgbFeatureVector {
    double values[10];
};

struct MethylXgbPrediction {
    double probability;
    bool somatic;
};

MethylXgbPrediction predictMethylXgb(MethylXgbVariantKind kind, const MethylXgbFeatureVector &features, double threshold);

#endif
