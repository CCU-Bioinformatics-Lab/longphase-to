#ifndef SOMATIC_REFINEMENT_POLICY_H
#define SOMATIC_REFINEMENT_POLICY_H

namespace somatic_refinement {

constexpr double kMethylXgbMaxPurity = 0.7;
constexpr double kHighPurityThreshold = 0.9;

struct Plan {
    bool runMethylXgb;
    bool convertNonGermlineToSomatic;
    bool rerunPhasing;
};

inline Plan makeResolvedPlan(double purity, bool methylXgbEnabled) {
    const bool runMethylXgb = methylXgbEnabled && purity <= kMethylXgbMaxPurity;
    const bool convertNonGermlineToSomatic = purity > kHighPurityThreshold;
    return {runMethylXgb, convertNonGermlineToSomatic,
            runMethylXgb || convertNonGermlineToSomatic};
}

inline bool shouldCollectMethylCalls(double configuredPurity, bool methylXgbEnabled) {
    return methylXgbEnabled &&
           (configuredPurity < 0.0 || configuredPurity <= kMethylXgbMaxPurity);
}

}  // namespace somatic_refinement

#endif
