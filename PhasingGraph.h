#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBaseMap;
// rrr, rra, rar, raa, arr, ara, arr, aar, aaa
using ThreePointEdge = std::array<float, 2*2*2>;

enum VariantOriginType{
    GERMLINE = 0,
    SOMATIC = 1,
    Origin_UNDEFINED = -1
};

enum EdgeType{
    edge_rrr = 0,
    edge_rra = 1,
    edge_rar = 2,
    edge_raa = 3,
    edge_arr = 4,
    edge_ara = 5,
    edge_aar = 6,
    edge_aaa = 7
};

// "somatic" corresponds to variable `a`
// "confidence" corresponds to variable `b`
// "source haplotype(allele)" corresponds to variable `c`
// `a_b_c` is used for their combined representation
enum SomaticVote{
    VOTE_UNDEFINED = 0,
    MID_HIGH_SR = 1,
    RIGHT_HIGH_SR = 2,
    LEFT_HIGH_SR = 3,

    MID_HIGH_SA = 5,
    RIGHT_HIGH_SA = 6,
    LEFT_HIGH_SA = 7,

    DISAGREE = 8,
    MID_LOW = 9,
    LEFT_LOW = 10,
    RIGHT_LOW = 11,
};

struct SomaticPattern {
    std::array<size_t, 3> edges;
    SomaticVote vote;
};

constexpr std::array<SomaticPattern, 12> highSomaticPatterns = {{
    {{edge_rrr, edge_ara, edge_aaa}, MID_HIGH_SA},
    {{edge_rra, edge_arr, edge_aar}, MID_HIGH_SA},
    {{edge_arr, edge_rra, edge_raa}, MID_HIGH_SR},
    {{edge_ara, edge_rrr, edge_rar}, MID_HIGH_SR},

    {{edge_rrr, edge_raa, edge_aaa}, LEFT_HIGH_SA},
    {{edge_rra, edge_rar, edge_aar}, LEFT_HIGH_SA},
    {{edge_rar, edge_rra, edge_ara}, LEFT_HIGH_SR},
    {{edge_raa, edge_rrr, edge_arr}, LEFT_HIGH_SR},

    {{edge_rrr, edge_aar, edge_aaa}, RIGHT_HIGH_SR},
    {{edge_arr, edge_rar, edge_raa}, RIGHT_HIGH_SR},
    {{edge_rar, edge_arr, edge_ara}, RIGHT_HIGH_SA},
    {{edge_aar, edge_rrr, edge_rra}, RIGHT_HIGH_SA},
}};

constexpr std::array<SomaticPattern, 3> lowSomaticPatterns = {{
    {{edge_ara, edge_ara, edge_aaa}, MID_LOW},
    {{edge_raa, edge_raa, edge_aaa}, LEFT_LOW},
    {{edge_aar, edge_aar, edge_aaa}, RIGHT_LOW},
}};

class SubEdge{
    
    private:
        int readCount;
        // Edge information. The vector store next pos
        // < next position, read name >
        std::map<int, std::vector<std::string> > *refRead;
        std::map<int, std::vector<std::string> > *altRead;
        // sum of edge pair quality, pos1 quality + pos2 quality
        // < next position, quality sum >
        std::map<int, int> *refQuality;
        std::map<int, int> *altQuality;
        // < next position, read count >
        std::map<int, float> *refReadCount;
        std::map<int, float> *altReadCount;

    public:
        
        SubEdge();
        ~SubEdge();
        
        void destroy();
        
        void addSubEdge(Variant &currentNode, Variant &connectNode, std::string readName, int baseQuality, double edgeWeight, bool fakeRead);
        std::pair<float,float> BestPair(int targetPos);
        float getRefReadCount(int targetPos);
        float getAltReadCount(int targetPos);        
        
        std::vector<std::string> showEdge(std::string message);
        std::vector<std::pair<int,int>> getConnectPos();
 
        int getQuality(PosAllele targetPos);
        int getAvgQuality(PosAllele targetPos);

};

//use to store the voting info from the previous variants
struct VoteResult{
    int Pos;		//who votes
    float para;		//rr+aa
    float cross;	//ra+ar
    float weight;	//how much weight
    int hap;		//which haplotype 
    double ESR;		//similarity of para and cross

    VoteResult( int currPos, float weight ) ;
};

struct VariantInfo {
    VariantType type;
    // bool homozygous; // 0 heterozygous 1 homozygous
    VariantOriginType origin = Origin_UNDEFINED; // 0 germline 1 somatic -1 unknown
    Allele sourceHaplotype = Allele_UNDEFINED; // REF_ALLELE ALT_ALLELE
    std::map<int, Allele> assignHaplotype;
};

struct VariantEdge{
    int currPos;
    SubEdge* alt;
    SubEdge* ref;
    std::map<int, std::map<int, ThreePointEdge>> *edgeCount;
    int refcnt ; // count the ref base amount
    int altcnt ; // count the alt base amount
    double vaf ; // count the vaf of the left snp
    int coverage ; // count the coverge on the snp
    
    VariantEdge(int currPos);
    // node pair 
    std::pair<PosAllele,PosAllele> findBestEdgePair(std::map<int, VariantInfo>::iterator currNodeIter, std::map<int, VariantInfo>::iterator nextNodeIter, bool isONT, double diffRatioThreshold, VoteResult &vote, bool debug);
    // number of read of two node. AA and AB combination
    std::pair<float,float> findNumberOfRead(int targetPos);
    // find the edge weight of two node
    // rrr, rra, rar, raa, arr, ara, aar, aaa
    ThreePointEdge findThreePointEdge(int pos, int nextPos);

    void addSegmentEdge(Variant &currentNode, Variant &connectNode, Variant &connectSecondNode);

    bool get_fakeSnp();
};


struct BlockRead{
    std::map<std::string,int> readVec;
    
    void recordRead(std::string readName);
};

struct EdgeResult{
    int rr;
    int ra;
    int ar;
    int aa;
};

class VairiantGraph{
    
    private:
        PhasingParameters *params;
        std::string *ref;
        std::string chr;
        std::vector<std::string> dotResult;
        std::vector<ReadVariant> *readVariant;
        
        // By default, a Map in C++ is sorted in increasing order based on its key.
        // position, edge
        std::map<int,VariantEdge*> *edgeList;
        // position, type
        std::map<int,VariantInfo> *variantPosType;

        // position phasing result
        PosPhasingResult *posPhasingResult;
        // store phased read and read's haplotype
        std::map<std::string,int> *readHpMap;

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection();

        void edgeConnectResult();

        int patternMining(ThreePointEdge threePointEdge);
        
        std::pair<float,float> Onelongcase( std::vector<VoteResult> vote ) ;

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params, std::string &chr);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> &in_readVariant);
        
        void phasingProcess(PosPhasingResult &posPhasingResult);

        void exportPhasingResult(PosPhasingResult &posPhasingResult);

        void somaticCalling(std::map<int, RefAlt>* variants);
        
        void writingDotFile(std::string dotPrefix);
        std::map<std::string,int>* getReadHP();
        int totalNode();

        void destroy();
        
};




#endif
