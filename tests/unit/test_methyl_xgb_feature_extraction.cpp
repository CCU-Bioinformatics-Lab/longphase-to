// Unit tests for methyl_xgb_feature_extraction::classifyAlleleAtAnchor.
//
// BAM SEQ (and pysam query_sequence) is stored in reference/forward orientation
// regardless of the read strand, and VCF ALT is forward-strand as well, so the
// allele classification must not depend on the read strand.
// See CCU-Bioinformatics-Lab/longphase-to#9.

#include "MethylXgbFeatureExtraction.h"

#include <iostream>
#include <string>
#include <vector>

using methyl_xgb_feature_extraction::AlleleObservation;
using methyl_xgb_feature_extraction::CigarOp;
using methyl_xgb_feature_extraction::classifyAlleleAtAnchor;

namespace {

int failures = 0;

const char *alleleName(int allele) {
    switch(allele) {
        case REF_ALLELE: return "REF";
        case ALT_ALLELE: return "ALT";
        default: return "UNDEFINED";
    }
}

// Read aligned at 0-based reference position 100; the variant anchor is 104
// (the last reference base before an insertion or deletion).
void expectAllele(const std::string &name,
                  const std::string &ref,
                  const std::string &alt,
                  const std::string &storedSequence,
                  const std::vector<CigarOp> &cigar,
                  int expected) {
    for(int reverse = 0; reverse <= 1; reverse++) {
        const AlleleObservation observed = classifyAlleleAtAnchor(
            100, 104, ref, alt, storedSequence, reverse == 1, cigar);
        const bool ok = observed.exactAllele == expected;
        std::cout << (ok ? "PASS " : "FAIL ") << name
                  << (reverse ? " [reverse]" : " [forward]")
                  << ": expected " << alleleName(expected)
                  << ", got " << alleleName(observed.exactAllele) << "\n";
        if(!ok) {
            failures++;
        }
    }
}

}  // namespace

int main() {
    const std::vector<CigarOp> match10 = {CigarOp('M', 10)};
    const std::vector<CigarOp> match12 = {CigarOp('M', 12)};
    const std::vector<CigarOp> insertion3 = {CigarOp('M', 5), CigarOp('I', 3), CigarOp('M', 5)};
    const std::vector<CigarOp> deletion2 = {CigarOp('M', 5), CigarOp('D', 2), CigarOp('M', 5)};

    expectAllele("SNV carrying ALT", "C", "T", "AAAATTTTTT", match10, ALT_ALLELE);
    expectAllele("SNV carrying REF", "C", "T", "AAAACTTTTT", match10, REF_ALLELE);

    // Inserted bases GTA stored in forward orientation on both strands.
    expectAllele("insertion carrying ALT", "C", "CGTA", "AAAACGTATTTTT", insertion3, ALT_ALLELE);
    // Inserted bases TAC are the reverse complement of GTA: must not match ALT
    // on either strand.
    expectAllele("insertion with other sequence", "C", "CGTA", "AAAACTACTTTTT", insertion3,
                 Allele_UNDEFINED);
    expectAllele("insertion absent", "C", "CGTA", "AAAACTTTTT", match10, REF_ALLELE);

    expectAllele("deletion carrying ALT", "CGG", "C", "AAAACTTTTT", deletion2, ALT_ALLELE);
    expectAllele("deletion absent", "CGG", "C", "AAAACGGTTTTT", match12, REF_ALLELE);

    if(failures != 0) {
        std::cout << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all checks passed\n";
    return 0;
}
