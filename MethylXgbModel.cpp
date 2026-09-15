#include "MethylXgbModel.h"

#include <cmath>

namespace {
struct MethylXgbNode {
    int feature;
    float threshold;
    int yes;
    int no;
    int missing;
    float leaf;
};

struct MethylXgbModel {
    const MethylXgbNode *nodes;
    const int *roots;
    int treeCount;
    double baseMargin;
};

#include "MethylXgbModelData.inc"

double sigmoid(double value) {
    if(value >= 0.0) {
        const double z = std::exp(-value);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(value);
    return z / (1.0 + z);
}

double predictRawMargin(const MethylXgbModel &model, const MethylXgbFeatureVector &features) {
    float values[10];
    for(int i = 0; i < 10; i++) {
        values[i] = static_cast<float>(features.values[i]);
    }

    double margin = model.baseMargin;
    for(int treeIdx = 0; treeIdx < model.treeCount; treeIdx++) {
        int nodeIdx = model.roots[treeIdx];
        while(true) {
            const MethylXgbNode &node = model.nodes[nodeIdx];
            if(node.feature < 0) {
                margin += node.leaf;
                break;
            }

            const float value = values[node.feature];
            if(std::isnan(value)) {
                nodeIdx = node.missing;
            }
            else if(value < node.threshold) {
                nodeIdx = node.yes;
            }
            else {
                nodeIdx = node.no;
            }
        }
    }
    return margin;
}

const MethylXgbModel &selectModel(MethylXgbVariantKind kind) {
    return kind == METHYL_XGB_KIND_SNV ? kSnvModel : kIndelModel;
}
}

MethylXgbPrediction predictMethylXgb(MethylXgbVariantKind kind, const MethylXgbFeatureVector &features, double threshold) {
    const MethylXgbModel &model = selectModel(kind);
    const double probability = sigmoid(predictRawMargin(model, features));

    MethylXgbPrediction prediction;
    prediction.probability = probability;
    prediction.somatic = probability >= threshold;
    return prediction;
}
