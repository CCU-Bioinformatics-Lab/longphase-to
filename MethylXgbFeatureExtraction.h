#ifndef METHYLXGBFEATUREEXTRACTION_H
#define METHYLXGBFEATUREEXTRACTION_H

#include "Util.h"

#include <string>
#include <vector>

namespace methyl_xgb_feature_extraction {

struct CigarOp {
    CigarOp(char inOp, int inLength) : op(inOp), length(inLength) {}

    char op;
    int length;
};

struct AlleleObservation {
    AlleleObservation(int inAllele = Allele_UNDEFINED,
                      bool inRawDetailEligible = false)
        : exactAllele(inAllele), rawDetailEligible(inRawDetailEligible) {}

    int exactAllele;
    bool rawDetailEligible;
};

// Reproduces the training builders' pileup-anchor allele classification.
// Coordinates are 0-based. Flags, MAPQ, and tumor-BAM eligibility are handled
// by the caller; rawDetailEligible only represents the legacy extractor's
// anchor-base screen.
AlleleObservation classifyAlleleAtAnchor(
    int alignmentStart0,
    int anchor0,
    const std::string &ref,
    const std::string &alt,
    const std::string &querySequenceAsStored,
    bool isReverse,
    const std::vector<CigarOp> &cigar);

// The training compact tables wrote each count ratio with Python-style .6f
// formatting and then read it back before computing site-level features.
double ratioAtTrainingPrecision(int numerator, int denominator);

}  // namespace methyl_xgb_feature_extraction

#endif
