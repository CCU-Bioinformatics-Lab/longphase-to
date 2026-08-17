#include "MethylXgbFeatureExtraction.h"

#include <cctype>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace methyl_xgb_feature_extraction {
namespace {

struct AnchorContext {
    AnchorContext() : mapped(false), queryPosition(0), followingIndel(0) {}

    bool mapped;
    std::size_t queryPosition;
    int followingIndel;
};

bool consumesReference(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N';
}

bool consumesQuery(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'S';
}

bool isAlignedBaseOp(char op) {
    return op == 'M' || op == '=' || op == 'X';
}

bool isSupportedCigarOp(char op) {
    return consumesReference(op) || consumesQuery(op) || op == 'H' || op == 'P';
}

AnchorContext findAnchorContext(int alignmentStart0,
                                int anchor0,
                                const std::vector<CigarOp> &cigar) {
    AnchorContext context;
    if(alignmentStart0 < 0 || anchor0 < alignmentStart0) {
        return context;
    }

    long long referencePosition = alignmentStart0;
    long long queryPosition = 0;

    for(std::size_t cigarIndex = 0; cigarIndex < cigar.size(); cigarIndex++) {
        const CigarOp &current = cigar[cigarIndex];
        if(current.length <= 0 || !isSupportedCigarOp(current.op)) {
            return AnchorContext();
        }

        const long long referenceEnd = referencePosition + current.length;
        if(isAlignedBaseOp(current.op) &&
           anchor0 >= referencePosition && anchor0 < referenceEnd) {
            context.mapped = true;
            context.queryPosition = static_cast<std::size_t>(
                queryPosition + (anchor0 - referencePosition));

            // pysam PileupRead.indel describes an I/D immediately after the
            // current pileup base. It is nonzero only at the final reference
            // base of this aligned CIGAR operation.
            if(anchor0 == referenceEnd - 1 && cigarIndex + 1 < cigar.size()) {
                const CigarOp &next = cigar[cigarIndex + 1];
                if(next.length <= 0 || !isSupportedCigarOp(next.op)) {
                    return AnchorContext();
                }
                if(next.op == 'I') {
                    context.followingIndel = next.length;
                }
                else if(next.op == 'D') {
                    context.followingIndel = -next.length;
                }
            }
            return context;
        }

        // A pileup anchor inside D/N has no query position.
        if((current.op == 'D' || current.op == 'N') &&
           anchor0 >= referencePosition && anchor0 < referenceEnd) {
            return AnchorContext();
        }

        if(consumesReference(current.op)) {
            referencePosition = referenceEnd;
        }
        if(consumesQuery(current.op)) {
            queryPosition += current.length;
        }
    }

    return context;
}

char uppercase(char base) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
}

char complement(char base) {
    switch(base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'n': return 'n';
        default: return base;
    }
}

std::string reverseComplement(const std::string &sequence) {
    std::string result;
    result.reserve(sequence.size());
    for(std::string::const_reverse_iterator iter = sequence.rbegin();
        iter != sequence.rend(); ++iter) {
        result.push_back(complement(*iter));
    }
    return result;
}

}  // namespace

AlleleObservation classifyAlleleAtAnchor(
    int alignmentStart0,
    int anchor0,
    const std::string &ref,
    const std::string &alt,
    const std::string &querySequenceAsStored,
    bool isReverse,
    const std::vector<CigarOp> &cigar) {
    if(ref.empty() || alt.empty()) {
        return AlleleObservation();
    }

    const AnchorContext context = findAnchorContext(alignmentStart0, anchor0, cigar);
    if(!context.mapped || context.queryPosition >= querySequenceAsStored.size()) {
        return AlleleObservation();
    }

    const char anchorBase = uppercase(querySequenceAsStored[context.queryPosition]);
    const bool rawDetailEligible =
        anchorBase == uppercase(ref[0]) || anchorBase == uppercase(alt[0]);

    if(ref.size() == 1 && alt.size() == 1) {
        if(anchorBase == uppercase(ref[0])) {
            return AlleleObservation(REF_ALLELE, rawDetailEligible);
        }
        if(anchorBase == uppercase(alt[0])) {
            return AlleleObservation(ALT_ALLELE, rawDetailEligible);
        }
        return AlleleObservation(Allele_UNDEFINED, rawDetailEligible);
    }

    if(alt.size() > ref.size()) {
        const std::size_t insertionLength = alt.size() - ref.size();
        if(context.followingIndel == static_cast<int>(insertionLength)) {
            const std::size_t insertionStart = context.queryPosition + 1;
            if(insertionStart > querySequenceAsStored.size() ||
               insertionLength > querySequenceAsStored.size() - insertionStart) {
                return AlleleObservation(Allele_UNDEFINED, rawDetailEligible);
            }

            std::string inserted = querySequenceAsStored.substr(
                insertionStart, insertionLength);
            if(isReverse) {
                inserted = reverseComplement(inserted);
            }
            if(inserted == alt.substr(ref.size())) {
                return AlleleObservation(ALT_ALLELE, rawDetailEligible);
            }
            return AlleleObservation(Allele_UNDEFINED, rawDetailEligible);
        }
        if(context.followingIndel == 0) {
            return AlleleObservation(REF_ALLELE, rawDetailEligible);
        }
        return AlleleObservation(Allele_UNDEFINED, rawDetailEligible);
    }

    if(ref.size() > alt.size()) {
        const std::size_t deletionLength = ref.size() - alt.size();
        if(context.followingIndel == -static_cast<int>(deletionLength)) {
            return AlleleObservation(ALT_ALLELE, rawDetailEligible);
        }
        if(context.followingIndel == 0) {
            return AlleleObservation(REF_ALLELE, rawDetailEligible);
        }
        return AlleleObservation(Allele_UNDEFINED, rawDetailEligible);
    }

    // Equal-length multi-base substitutions were not classified by the
    // deployment builders.
    return AlleleObservation(Allele_UNDEFINED, rawDetailEligible);
}

double ratioAtTrainingPrecision(int numerator, int denominator) {
    if(numerator < 0 || denominator <= 0 || numerator > denominator) {
        throw std::invalid_argument("invalid methylation ratio counts");
    }

    const double ratio = static_cast<double>(numerator) /
                         static_cast<double>(denominator);
    std::ostringstream formatted;
    formatted.imbue(std::locale::classic());
    formatted << std::fixed << std::setprecision(6) << ratio;

    double result = 0.0;
    std::istringstream parsed(formatted.str());
    parsed.imbue(std::locale::classic());
    parsed >> result;
    if(!parsed) {
        throw std::runtime_error("failed to parse six-decimal methylation ratio");
    }
    return result;
}

}  // namespace methyl_xgb_feature_extraction
