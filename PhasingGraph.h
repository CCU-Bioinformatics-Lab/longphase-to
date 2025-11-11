#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBaseMap;
// rrr, rra, rar, raa, arr, ara, arr, aar, aaa
using ThreePointEdge = std::array<float, 2*2*2>;

enum EdgeType{
    EDGE_RRR = 0,
    EDGE_RRA = 1,
    EDGE_RAR = 2,
    EDGE_RAA = 3,
    EDGE_ARR = 4,
    EDGE_ARA = 5,
    EDGE_AAR = 6,
    EDGE_AAA = 7
};

enum SomaticVote{
    VOTE_UNDEFINED = 0,
    DISAGREE = 1,

    MID_HIGH_SA_PAR = 2,
    MID_HIGH_SR_PAR = 3,
    MID_HIGH_SR_CR = 4,
    MID_HIGH_SA_CR = 5,

    RIGHT_HIGH_SR_PAR = 6,
    RIGHT_HIGH_SA_PAR = 7,
    RIGHT_HIGH_SA_CR = 8,
    RIGHT_HIGH_SR_CR = 9,

    LEFT_HIGH_SA_PAR = 10,
    LEFT_HIGH_SR_PAR = 11,
    LEFT_HIGH_SR_CR = 12,
    LEFT_HIGH_SA_CR = 13,

    MID_LOW_1 = 14,
    MID_LOW_2 = 15,
    MID_LOW_3 = 16,

    LEFT_LOW_1 = 17,
    LEFT_LOW_2 = 18,
    LEFT_LOW_3 = 19,

    RIGHT_LOW_1 = 20,
    RIGHT_LOW_2 = 21,
    RIGHT_LOW_3 = 22,
};

struct SomaticPattern {
    std::array<size_t, 3> edges;
    SomaticVote vote;
};

constexpr std::array<SomaticPattern, 12> highSomaticPatterns = {{
    {{EDGE_RRR, EDGE_ARA, EDGE_AAA}, MID_HIGH_SA_PAR},
    {{EDGE_RRA, EDGE_ARR, EDGE_AAR}, MID_HIGH_SA_CR},
    {{EDGE_ARR, EDGE_RRA, EDGE_RAA}, MID_HIGH_SR_CR},
    {{EDGE_ARA, EDGE_RRR, EDGE_RAR}, MID_HIGH_SR_PAR},

    {{EDGE_RRR, EDGE_RAA, EDGE_AAA}, LEFT_HIGH_SA_PAR},
    {{EDGE_RRA, EDGE_RAR, EDGE_AAR}, LEFT_HIGH_SA_CR},
    {{EDGE_RAR, EDGE_RRA, EDGE_ARA}, LEFT_HIGH_SR_CR},
    {{EDGE_RAA, EDGE_RRR, EDGE_ARR}, LEFT_HIGH_SR_PAR},

    {{EDGE_RRR, EDGE_AAR, EDGE_AAA}, RIGHT_HIGH_SR_PAR},
    {{EDGE_ARR, EDGE_RAR, EDGE_RAA}, RIGHT_HIGH_SR_CR},
    {{EDGE_RAR, EDGE_ARR, EDGE_ARA}, RIGHT_HIGH_SA_CR},
    {{EDGE_AAR, EDGE_RRR, EDGE_RRA}, RIGHT_HIGH_SA_PAR},
}};

constexpr std::array<SomaticPattern, 9> lowSomaticPatterns = {{
    {{EDGE_ARA, EDGE_AAA, EDGE_RRR}, MID_LOW_1},
    {{EDGE_ARA, EDGE_AAA, EDGE_RRA}, MID_LOW_2},
    {{EDGE_ARA, EDGE_AAA, EDGE_ARR}, MID_LOW_3},

    {{EDGE_RAA, EDGE_AAA, EDGE_RRR}, LEFT_LOW_1},
    {{EDGE_RAA, EDGE_AAA, EDGE_RRA}, LEFT_LOW_2},
    {{EDGE_RAA, EDGE_AAA, EDGE_RAR}, LEFT_LOW_3},

    {{EDGE_AAR, EDGE_AAA, EDGE_RRR}, RIGHT_LOW_1},
    {{EDGE_AAR, EDGE_AAA, EDGE_RAR}, RIGHT_LOW_2},
    {{EDGE_AAR, EDGE_AAA, EDGE_ARR}, RIGHT_LOW_3},
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
    bool homozygous; // 0 heterozygous 1 homozygous
    VariantOriginType origin = ORIGIN_UNDEFINED; // 0 pon 1 somatic 2 germline -1 unknown
};

struct VariantEdge{
    int currPos;
    SubEdge* alt;
    SubEdge* ref;
    std::map<int, std::map<int, ThreePointEdge>> *edgeCount;
    int refcnt ; // count the ref base amount
    int altcnt ; // count the alt base amount
    double vaf = 0; // count the vaf of the left snp
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

struct VariantBases{
    int targetCount;
    std::map<int, int> offsetDiffRefCount; //offset, diff ref count
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

        // Map to store distribution statistics
        std::map<double, int> *ploidyRatioMap;

        // Process read variants to determine haplotype and update allele counts
        void processReadVariants(std::map<int,std::map<int,std::map<double,double>>> *hpAlleleCountMap);

        // Calculate reference allele distribution
        void calculatePloidyRatio(double hp1Ref, double hp2Ref, std::map<double, int> *ploidyRatioMap);

        // Reassign allele results based on haplotype counts
        void reassignAlleleResult(std::map<int,std::map<int,std::map<double,double>>> *hpAlleleCountMap, std::map<double, int> *ploidyRatioMap);

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection(std::map<double, int> *ploidyRatioMap);

        void edgeConnectResult(std::vector<LOHSegment> &LOHSegments);

        SomaticVote patternMining(ThreePointEdge threePointEdge, bool leftPosIndel, bool middlePosIndel, bool rightPosIndel);
        
        std::pair<float,float> Onelongcase( std::vector<VoteResult> vote ) ;

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params, std::string &chr);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> *in_readVariant);
        
        void phasingProcess(PosPhasingResult &posPhasingResult, std::vector<LOHSegment> &LOHSegments, std::map<double, int> *ploidyRatioMap);

        void somaticCalling(std::map<int, RefAlt>* variants);

        void tagSomatic(std::map<int, RefAlt>* variants);

        void convertNonGermlineToSomatic();

        void exportPhasingResult(PosPhasingResult &posPhasingResult, std::vector<LOHSegment> &LOHSegments);
        
        void writingDotFile(std::string dotPrefix);
        std::map<std::string,int>* getReadHP();
        int totalNode();

        void destroy();
        
};

class Roller {
    public:
        Roller();
        Roller(size_t windowSize);
        Roller(size_t windowSize, const std::vector<double>& input);
        
        Roller& setValues(const std::vector<double>& input);
        Roller& smooth();
        Roller& forward();
        Roller& reverse();
        Roller& downSampling(int distance);
        Roller& directionalDifference(size_t distance = 10);
        Roller& opposite();
        Roller& reverseNetValues();
        // Get the final result
        std::vector<double> getResult() const;
        
        std::vector<double> smooth(const std::vector<double>& values);
        std::vector<double> forward(const std::vector<double>& values);
        std::vector<double> reverse(const std::vector<double>& values);
        std::vector<double> downSampling(const std::vector<double>& values, int distance);
        std::vector<double> directionalDifference(const std::vector<double>& values, size_t distance = 10);
        std::vector<double> opposite(const std::vector<double>& values);
        std::vector<double> reverseNetValues(const std::vector<double>& values);

        static std::vector<double> backValues(const std::vector<size_t>& idxs, const std::vector<double>& values);
        static std::vector<int> backPosition(const std::vector<size_t>& idxs, const std::vector<int>& values);
        static void replaceValue(ClipCount &clipCount, std::vector<int>& keys, std::vector<size_t>& idxs, std::vector<double>& replaceValues);

    private:
        size_t window;
        std::vector<double> values;
        void updateWindow(double& rollingSum, std::deque<double>& windowValues, double newValue);
};

class PointFinder {
    public:
        PointFinder(size_t peakDistance, int calibrationThreshold = 5, size_t calibrationDistance = 0);
        
        std::vector<size_t> findAllPeaks(const std::vector<double>& values, double peakThreshold = 0);
        std::vector<size_t> findPeaks(const std::vector<double>& values, double peakThreshold = 0);
        std::vector<size_t> getForwardGentle(const std::vector<double>& values, std::vector<size_t> peaks);
        std::vector<size_t> getReverseGentle(const std::vector<double>& values, std::vector<size_t> peaks);

    private:
        size_t peakDistance;
        int calibrationThreshold;
        size_t calibrationDistance;

        size_t start(size_t target, size_t distance, size_t distanceDivisor = 1);
        size_t end(size_t target, size_t distance, size_t size, size_t distanceDivisor = 1);
};

class Clip{
    private:
        std::string chr;
        std::vector<int> *largeGenomicEventInterval;
        std::vector<std::pair<int, int>> *smallGenomicEventRegion;
        void amplifyGentleClip(ClipCount &clipCount);
        void detectInterval(ClipCount &clipCount);
    public:
        Clip(std::string &chr);
        ~Clip();
        void detectGenomicEventInterval(ClipCount &clipCount, std::vector<int> &largeGenomicEventInterval, std::vector<std::pair<int, int>> &smallGenomicEventRegion);
        void detectLOHRegion(SnpParser &snpMap, std::vector<LOHSegment> &LOHSegments);
};

class PurityCalculator {
    private:
        static std::map<double, int> mergeDistributionMap(const std::map<std::string, std::map<double, int>>& data);
        static int getTotalCount(const std::map<double, int>& data);
        static double findQuartile(const std::map<double, int>& data, double targetPos);
        static double getLOHRatio(std::map<std::string, ChrInfo> &chrInfo, std::map<std::string, int> &chrLength);
    public:
        static double getPurity(std::map<std::string, std::map<double, int>> &inChrDistributionMap, std::string &output_root_path, Caller caller, std::map<std::string, ChrInfo> &chrInfo, std::map<std::string, int> &chrLength);
};

#endif
