#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <vector>
#include <omp.h>
#include <deque>

enum Haplotype {
    HAPLOTYPE_UNDEFINED = -1,
    HAPLOTYPE1 = 0,
    HAPLOTYPE2 = 1
};

enum Allele {
    Allele_UNDEFINED = -1,
    REF_ALLELE = 0,
    ALT_ALLELE = 1
};

enum VariantType {
    // A VariantType value > 0 represents the quality of the read variant base and indicates that its type is SNP = -2.
    VARIANT_UNDEFINED = -1, // Undefined variant type
    SNP = -2, // Single Nucleotide Polymorphism
    INDEL = -3, // Insertion/Deletion
    SV = -4, // Structure Variation
    MOD_FORWARD_STRAND = -5, // Forward Strand Modification
    MOD_REVERSE_STRAND = -6, // Reverse Strand Modification
    DANGER_INDEL = -7, // Danger Insertion/Deletion
};

enum VariantGenotype {
    HOM = 1, // homozygous
    HET = 0, // heterozygous
    GENOTYPE_UNDEFINED = -1, // undefined genotype
};

enum ClipFrontBack {
    FRONT = 0,
    BACK = 1
};

constexpr size_t PHASING_RESULT_SIZE = 3;
struct PhasingResult {
    Haplotype refHaplotype = HAPLOTYPE_UNDEFINED;
    std::vector<int> phaseSet;
    VariantType type = VARIANT_UNDEFINED;
    bool somatic = false;
    std::vector<std::string> genotype;
    PhasingResult() = default;
    PhasingResult(Haplotype inRefHaplotype, int inPhaseSet, VariantType inType)
        : refHaplotype(inRefHaplotype), type(inType){
            phaseSet.push_back(inPhaseSet);
    }
    PhasingResult(Haplotype inRefHaplotype, int inPhaseSet, VariantType inType, bool inSomatic)
        : refHaplotype(inRefHaplotype), type(inType), somatic(inSomatic) {
            phaseSet.push_back(inPhaseSet);
    }
    PhasingResult(int inPhaseSet, VariantType inType, std::string inGenotype)
        : type(inType){
            phaseSet.push_back(inPhaseSet);
            genotype.push_back(inGenotype);
    }
};
typedef std::map<int, PhasingResult> PosPhasingResult;
// pos<read start|read end ,count >
typedef std::map<int, std::array<int, 2>> ClipCount;
struct LOHSegment{
    int start;
    int end;
    Allele startAllele;
    Allele endAllele;
    LOHSegment(int inStart, int inEnd): 
        start(inStart), end(inEnd), startAllele(Allele_UNDEFINED), endAllele(Allele_UNDEFINED){}
};
struct ChrInfo{
    PosPhasingResult posPhasingResult = PosPhasingResult();
    ClipCount clipCount = ClipCount();
    std::vector<int> largeGenomicEventInterval = std::vector<int>();
    std::vector<std::pair<int, int>> smallGenomicEventRegion = std::vector<std::pair<int, int>>();
    std::vector<LOHSegment> LOHSegments = std::vector<LOHSegment>();
};
typedef std::map<std::string, std::map<int, PhasingResult>> ChrPhasingResult;


// use for parsing
struct Variant{
    Variant(int position, int allele, int quality, bool homozygous = false):
    position(position), 
    allele(allele), 
    quality(quality),
    homozygous(homozygous){};
    
    int position;
    int allele;
    int quality;
    bool underHomopolymer;
    bool homozygous;
};

struct ReadVariant{
    // init function
    ReadVariant(): read_name(""), 
            mapping_quality(0), 
            source_id(""), 
            sample_id(""), 
            reference_start(0), 
            BX_tag(""){}
            
    std::string read_name;
    int mapping_quality;
    std::string source_id;
    std::string sample_id;
    int reference_start;
    std::string BX_tag;
    bool is_reverse;
    bool fakeRead;
    
    std::vector<Variant> variantVec;
    
    void sort();
};

struct less_than_key
{
    inline bool operator() (const Variant& v1, const Variant& v2)
    {
        return (v1.position < v2.position);
    }
};


std::string getTargetString(std::string line, std::string start_sign, std::string end_sign);

int homopolymerLength(int snp_pos, const std::string &ref_string);




#endif