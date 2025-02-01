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

enum Haplotype {
    HAPLOTYPE_UNDEFINED = -1,
    HAPLOTYPE1 = 0,
    HAPLOTYPE2 = 1
};

struct PhasingResult {
    int refHaplotype;
    int phaseSet;
    int type;

    PhasingResult(): refHaplotype(-1), phaseSet(-1), type(-1) {}
    PhasingResult(int inRefHaplotype, int inPhaseSet, int inType)
        : refHaplotype(inRefHaplotype), phaseSet(inPhaseSet), type(inType) {}
};
typedef std::map<int,PhasingResult> PosPhasingResult;
typedef std::map<std::string, PosPhasingResult> ChrPhasingResult;

// use for parsing
struct Variant{
    Variant(int position, int allele, int quality):
    position(position), 
    allele(allele), 
    quality(quality){};
    
    int position;
    int allele;
    int quality;
    bool underHomopolymer;
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