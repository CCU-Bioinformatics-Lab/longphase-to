#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>

constexpr size_t HAPLOTYPE_SIZE = 5;

struct HaplotagParameters
{
    int numThreads;
    int qualityThreshold;
    
    double percentageThreshold;
    
    std::string snpFile;
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string fastaFile;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    std::string outputFormat;
    
    bool tagSupplementary;
    bool writeReadLog;
};

struct AlleleHaplotype
{
    std::string ref;
    std::string alt;
    Haplotype refHaplotype = HAPLOTYPE_UNDEFINED;
    Haplotype altHaplotype = HAPLOTYPE_UNDEFINED;
};

class HaplotagProcess
{
    void variantParser(std::string variantFile);
    void compressParser(std::string &variantFile);
    void unCompressParser(std::string &variantFile);
    void parserProcess(std::string &input);
    
    void tagRead(HaplotagParameters &params);
    
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;
    
    // chr, variant position (0-base), allele haplotype
    std::map<std::string, std::map<int, AlleleHaplotype> > chrVariant;
    // chr, variant position (0-base), phased set
    std::map<std::string, int > psIndex;
    std::map<std::string, std::map<int, int> > chrVariantPS;
    
    std::map<int, AlleleHaplotype> currentChrVariants;
    std::map<int, AlleleHaplotype>::iterator firstVariantIter;
    // The number of SVs occurring on different haplotypes in a read
    std::map<std::string, std::map<int, int> > readSVHapCount;

    void initFlag(bam1_t *aln, std::string flag);
    
    std::string judgeHaplotype(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, const std::string &ref_string);
    void countSNPHaplotype(std::string &base, AlleleHaplotype &haplotypeBase, std::array<int, HAPLOTYPE_SIZE> &countMap, int & HP);
    void countINDELHaplotype(bool isRef, AlleleHaplotype &haplotypeBase, std::array<int, HAPLOTYPE_SIZE> &countMap, int & HP);
    void getVote(std::array<int, HAPLOTYPE_SIZE> &countMap, double &min, double &max, std::string &hpResult);
    
    int totalAlignment;
    int totalSupplementary;
    int totalSecondary;
    int totalUnmapped;
    int totalTagCount;
    int totalUnTagCount;
    
    std::time_t processBegin;
    bool integerPS;
    bool parseSnpFile;
    bool parseSVFile;
    bool parseMODFile;
    
    public:
        HaplotagProcess(HaplotagParameters params);
        ~HaplotagProcess();

};


#endif