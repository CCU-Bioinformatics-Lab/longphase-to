#include "Phasing.h"
#include "PhasingProcess.h"
#include "SomaticRefinementPolicy.h"
#include "Util.h"
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <limits>


#define SUBPROGRAM "phase"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"   --help                                 display this help and exit.\n\n"
"require arguments:\n"
"   -s, --snp-file=NAME                    input SNP vcf file, or a VCF file with both SNPs and indels.\n"
"   -b, --bam-file=NAME                    input bam file.\n"
"   -r, --reference=NAME                   reference fasta.\n"
"   -c, --caller=NAME                      variant caller name.\n"
"                                          options: clairs_to_ss, clairs_to_ssrs, deepsomatic_to\n\n"

"optional arguments:\n"
"   --sv-file=NAME                         input SV vcf file.\n"
"   --mod-file=NAME                        input modified vcf file.(produce by longphase modcall)\n"
"   -t, --threads=Num                      number of thread. default:1\n"
"   -o, --out-prefix=NAME                  prefix of phasing result. default: result\n"
"   --indels                               phase small indel. default: False\n"
"   --dot                                  each contig/chromosome will generate dot file.\n"
"   --loh                                  output LOH results. default: False\n\n"

"somatic arguments:\n"
"   --purity=[0~1]                         set sample purity directly; if omitted, estimate automatically.\n"
"   --disable-calling                      disable longphase calling mode. default: False\n"
"   --disable-pon-tag                      disable reading the VCF FILTER field PON to determine germline variants. default: False\n"
"   --disable-refine-somatic               do not modify VCF FILTER based on somatic refinement. default: False\n"
"   --pon-file=NAME                        input PON VCF file. determines germline variants using position-based matching.\n"
"                                          input format: A.vcf,B.vcf\n"
"   --strict-pon-file=NAME                 input PON VCF file. determines germline variants using both position and ALT allele matching.\n"
"                                          input format: A.vcf,B.vcf\n"
"   --somaticConnectAdjacent=Num           connect adjacent N SNPs. default:6\n\n"

"methylation XGBoost somatic refinement arguments:\n"
"   --methyl-xgb                           enable refinement when purity <=0.7. default: True\n"
"   --disable-methyl-xgb                   disable methylation feature XGBoost somatic refinement.\n"
"   --methyl-xgb-snv-threshold=[0~1]       SNV somatic probability threshold. default:0.44\n"
"   --methyl-xgb-indel-threshold=[0~1]     indel somatic probability threshold. default:0.17\n"
"   --methyl-window=Num                    variant-centered methylation window radius. default:2000\n"
"   --meth-high=[0~1]                      high methylation probability threshold. default:0.8\n"
"   --meth-low=[0~1]                       low methylation probability threshold. default:0.2\n\n"

"parse alignment arguments:\n"
"   -q, --mappingQuality=Num               filter alignment if mapping quality is lower than threshold. default:1\n"
"   -x, --mismatchRate=Num                 mark reads as false if mismatchRate of them are higher than threshold. default:3\n\n"

"phasing graph arguments:\n"
"   -p, --baseQuality=[0~90]               change edge's weight to --edgeWeight if base quality is lower than the threshold. default:12\n"
"   -e, --edgeWeight=[0~1]                 if one of the bases connected by the edge has a quality lower than --baseQuality\n"
"                                          its weight is reduced from the normal 1. default:0.1\n"
"   -a, --connectAdjacent=Num              connect adjacent N SNPs. default:35\n"
"   -d, --distance=Num                     phasing two variant if distance less than threshold. default:300000\n"
"   -1, --edgeThreshold=[0~1]              give up SNP-SNP phasing pair if the number of reads of the \n"
"                                          two combinations are similar. default:0.7\n"
"   -L, --overlapThreshold=[0~1]           filtering different alignments of the same read if there is overlap. default:0.2 \n\n"


"haplotag read correction arguments:\n"
"   -m, --readConfidence=[0.5~1]           The confidence of a read being assigned to any haplotype. default:0.65\n"
"   -n, --snpConfidence=[0.5~1]            The confidence of assigning two alleles of a SNP to different haplotypes. default:0.75\n\n";

static const char* shortopts = "s:b:o:t:r:d:1:a:q:x:p:e:n:m:L:c:";

enum { OPT_HELP = 1 , DOT_FILE, SV_FILE, MOD_FILE, IS_ONT, IS_PB, PHASE_INDEL, VERSION, PON_FILE, STRICT_PON_FILE, SOMATIC_CONNECT_ADJACENT, OUTPUT_LOH, OUTPUT_SGE, OUTPUT_LGE, OUTPUT_GE, DISABLE_PON_TAG, DISABLE_CALLING, DISABLE_REFINE_SOMATIC, OPT_PURITY, METHYL_XGB, DISABLE_METHYL_XGB, METHYL_XGB_SNV_THRESHOLD, METHYL_XGB_INDEL_THRESHOLD, METHYL_WINDOW, METH_HIGH, METH_LOW};

static const struct option longopts[] = {
    { "help",                 no_argument,        NULL, OPT_HELP },
    { "dot",                  no_argument,        NULL, DOT_FILE },  
    { "ont",                  no_argument,        NULL, IS_ONT }, 
    { "pb",                   no_argument,        NULL, IS_PB }, 
    { "version",              no_argument,        NULL, VERSION }, 
    { "indels",               no_argument,        NULL, PHASE_INDEL },   
    { "loh",                  no_argument,        NULL, OUTPUT_LOH },
    { "sge",                  no_argument,        NULL, OUTPUT_SGE },
    { "lge",                  no_argument,        NULL, OUTPUT_LGE },
    { "ge",                   no_argument,        NULL, OUTPUT_GE },
    { "sv-file",              required_argument,  NULL, SV_FILE },  
    { "mod-file",             required_argument,  NULL, MOD_FILE },
    { "pon-file",             required_argument,  NULL, PON_FILE },
    { "strict-pon-file",      required_argument,  NULL, STRICT_PON_FILE },
    { "somaticConnectAdjacent", required_argument,  NULL, SOMATIC_CONNECT_ADJACENT },
    { "disable-pon-tag",      no_argument,        NULL, DISABLE_PON_TAG },
    { "disable-calling",      no_argument,        NULL, DISABLE_CALLING },
    { "disable-refine-somatic", no_argument,        NULL, DISABLE_REFINE_SOMATIC },
    { "methyl-xgb",           no_argument,        NULL, METHYL_XGB },
    { "disable-methyl-xgb",   no_argument,        NULL, DISABLE_METHYL_XGB },
    { "methyl-xgb-snv-threshold", required_argument,  NULL, METHYL_XGB_SNV_THRESHOLD },
    { "methyl-xgb-indel-threshold", required_argument,  NULL, METHYL_XGB_INDEL_THRESHOLD },
    { "methyl-window",        required_argument,  NULL, METHYL_WINDOW },
    { "meth-high",            required_argument,  NULL, METH_HIGH },
    { "meth-low",             required_argument,  NULL, METH_LOW },
    { "reference",            required_argument,  NULL, 'r' },
    { "snp-file",             required_argument,  NULL, 's' },
    { "bam-file",             required_argument,  NULL, 'b' },
    { "out-prefix",           required_argument,  NULL, 'o' },
    { "threads",              required_argument,  NULL, 't' },
    { "distance",             required_argument,  NULL, 'd' },
    { "edgeThreshold",        required_argument,  NULL, '1' },
    { "connectAdjacent",      required_argument,  NULL, 'a' },
    { "mappingQuality",       required_argument,  NULL, 'q' },
    { "mismatchRate",       required_argument,  NULL, 'x' },
    { "baseQuality",       required_argument,  NULL, 'p' },
    { "edgeWeight",       required_argument,  NULL, 'e' },
    { "snpConfidence",        required_argument,  NULL, 'n' },
    { "readConfidence",       required_argument,  NULL, 'm' },
    { "overlapThreshold",     required_argument,  NULL, 'L' },
    { "caller",               required_argument,  NULL, 'c' },
    { "purity",               required_argument,  NULL, OPT_PURITY },
    { NULL, 0, NULL, 0 }
};

namespace {
bool hasOnlyTrailingWhitespace(const char *text) {
    while(text != nullptr && *text != '\0') {
        if(!std::isspace(static_cast<unsigned char>(*text))) {
            return false;
        }
        text++;
    }
    return true;
}

bool parseFiniteDouble(const char *text, double &value) {
    if(text == nullptr || *text == '\0') {
        return false;
    }

    errno = 0;
    char *end = nullptr;
    const double parsed = std::strtod(text, &end);
    if(end == text || errno == ERANGE || !hasOnlyTrailingWhitespace(end) || !std::isfinite(parsed)) {
        return false;
    }

    value = parsed;
    return true;
}

bool parseFiniteFloat(const char *text, float &value) {
    double parsed = 0.0;
    if(!parseFiniteDouble(text, parsed) ||
       parsed < -std::numeric_limits<float>::max() ||
       parsed > std::numeric_limits<float>::max()) {
        return false;
    }

    const float converted = static_cast<float>(parsed);
    if(!std::isfinite(converted)) {
        return false;
    }

    value = converted;
    return true;
}

bool parseInteger(const char *text, int &value) {
    if(text == nullptr || *text == '\0') {
        return false;
    }

    errno = 0;
    char *end = nullptr;
    const long parsed = std::strtol(text, &end, 10);
    if(end == text || errno == ERANGE || !hasOnlyTrailingWhitespace(end) ||
       parsed < std::numeric_limits<int>::min() ||
       parsed > std::numeric_limits<int>::max()) {
        return false;
    }

    value = static_cast<int>(parsed);
    return true;
}

const char *numericOptionValue(const char *text) {
    return text == nullptr ? "" : text;
}
}

namespace opt
{
    static int numThreads = 1;
    static int distance = 300000;
    static std::string snpFile="";
    static std::string svFile="";
    static std::string modFile="";
    static std::string ponFile="";
    static std::string strictPonFile="";
    static std::vector<std::string> bamFile;
    static std::string fastaFile="";
    static std::string resultPrefix="result";
    static std::string callerStr="";
    static Caller caller = CALLER_UNDEFINED;
    static bool generateDot=false;
    static bool phaseIndel=false;
    static bool disablePonTag=false;
    static bool disableCalling=false;
    static bool disableRefineSomatic=false;
    static bool enableMethylXgb=true;
    static bool methylXgbExplicitEnable=false;
    static bool methylXgbExplicitDisable=false;
    static double methylXgbSnvThreshold=METHYL_XGB_DEFAULT_SNV_THRESHOLD;
    static double methylXgbIndelThreshold=METHYL_XGB_DEFAULT_INDEL_THRESHOLD;
    static int methylXgbWindow=2000;
    static float methylXgbMethHigh=0.8f;
    static float methylXgbMethLow=0.2f;
    
    static int connectAdjacent = 35;
    static int mappingQuality = 1;
    static double mismatchRate = 3;
    
    static int baseQuality = 12;
    static double edgeWeight = 0.1 ;
    
    static double snpConfidence  = 0.75;
    static double readConfidence = 0.65;
    
    static double edgeThreshold = 0.7;

    static double overlapThreshold = 0.2;

    static std::string command;

    static int somaticConnectAdjacent = 6;

    static bool outputLOH = false;
    static bool outputSGE = false;
    static bool outputLGE = false;
    static bool outputGE = false;
    static double purity = -1.0; // user-provided purity; negative means unset
}

void PhasingOptions(int argc, char** argv)
{
    optind=1;    //reset getopt

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
        case 's': arg >> opt::snpFile; break;
        case 't': arg >> opt::numThreads; break;
        case 'o': arg >> opt::resultPrefix; break;
        case 'r': arg >> opt::fastaFile; break;  
        case 'd': arg >> opt::distance; break;  
        case '1': arg >> opt::edgeThreshold; break; 
        case 'a': arg >> opt::connectAdjacent; break;
        case 'q': arg >> opt::mappingQuality; break;
        case 'x': arg >> opt::mismatchRate; break;
        case 'p': arg >> opt::baseQuality; break;
        case 'e': arg >> opt::edgeWeight; break;
        case 'n': arg >> opt::snpConfidence; break;
        case 'm': arg >> opt::readConfidence; break;
        case 'L': arg >> opt::overlapThreshold; break;
        case 'c': arg >> opt::callerStr; break;
        case 'b': {
            std::string bamFile;
            arg >> bamFile;
            opt::bamFile.push_back(bamFile); break;
        }
        case SV_FILE:  arg >> opt::svFile; break; 
        case MOD_FILE: arg >> opt::modFile; break; 
        case PHASE_INDEL: opt::phaseIndel=true; break; 
        case DOT_FILE: opt::generateDot=true; break;
        case PON_FILE: arg >> opt::ponFile; break;
        case STRICT_PON_FILE: arg >> opt::strictPonFile; break;
        case SOMATIC_CONNECT_ADJACENT: arg >> opt::somaticConnectAdjacent; break;
        case OUTPUT_LOH: opt::outputLOH=true; break;
        case OUTPUT_SGE: opt::outputSGE=true; break;
        case OUTPUT_LGE: opt::outputLGE=true; break;
        case OUTPUT_GE: opt::outputGE=true; break;
        case DISABLE_PON_TAG: opt::disablePonTag=true; break;
        case DISABLE_CALLING: opt::disableCalling=true; break;
        case DISABLE_REFINE_SOMATIC: opt::disableRefineSomatic=true; break;
        case OPT_PURITY:
            if(!parseFiniteDouble(optarg, opt::purity) || opt::purity < 0.0 || opt::purity > 1.0) {
                std::cerr << SUBPROGRAM " invalid purity. value: "
                          << numericOptionValue(optarg)
                          << "\n please check --purity=[0~1]\n";
                die = true;
            }
            break;
        case METHYL_XGB: opt::enableMethylXgb=true; opt::methylXgbExplicitEnable=true; break;
        case DISABLE_METHYL_XGB: opt::enableMethylXgb=false; opt::methylXgbExplicitDisable=true; break;
        case METHYL_XGB_SNV_THRESHOLD:
            if(!parseFiniteDouble(optarg, opt::methylXgbSnvThreshold)) {
                std::cerr << SUBPROGRAM " invalid methyl-xgb-snv-threshold. value: "
                          << numericOptionValue(optarg)
                          << "\n please check --methyl-xgb-snv-threshold=[0~1]\n";
                die = true;
            }
            break;
        case METHYL_XGB_INDEL_THRESHOLD:
            if(!parseFiniteDouble(optarg, opt::methylXgbIndelThreshold)) {
                std::cerr << SUBPROGRAM " invalid methyl-xgb-indel-threshold. value: "
                          << numericOptionValue(optarg)
                          << "\n please check --methyl-xgb-indel-threshold=[0~1]\n";
                die = true;
            }
            break;
        case METHYL_WINDOW:
            if(!parseInteger(optarg, opt::methylXgbWindow)) {
                std::cerr << SUBPROGRAM " invalid methyl-window. value: "
                          << numericOptionValue(optarg)
                          << "\n please check --methyl-window=Num\n";
                die = true;
            }
            break;
        case METH_HIGH:
            if(!parseFiniteFloat(optarg, opt::methylXgbMethHigh)) {
                std::cerr << SUBPROGRAM " invalid meth-high. value: "
                          << numericOptionValue(optarg)
                          << "\n please check --meth-high=[0~1]\n";
                die = true;
            }
            break;
        case METH_LOW:
            if(!parseFiniteFloat(optarg, opt::methylXgbMethLow)) {
                std::cerr << SUBPROGRAM " invalid meth-low. value: "
                          << numericOptionValue(optarg)
                          << "\n please check --meth-low=[0~1]\n";
                die = true;
            }
            break;
        case OPT_HELP:
            std::cout << CORRECT_USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        }
    }
    
    for(int i = 0; i < argc; ++i){
        opt::command.append(argv[i]);
        opt::command.append(" ");
    }

    if (argc - optind < 0 )
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    
    if( opt::snpFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::snpFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing SNP file.\n";
        die = true;
    }
    
    if( opt::fastaFile != "")
    {
        std::ifstream openFile( opt::fastaFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::fastaFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing reference.\n";
        die = true;
    }

    if(opt::bamFile.empty()){
        std::cerr << SUBPROGRAM ": missing BAM file.\n";
        die = true;
    }

    if(opt::bamFile.size() > 1 &&
       somatic_refinement::shouldCollectMethylCalls(opt::purity, opt::enableMethylXgb)){
        std::cerr << SUBPROGRAM
                  << ": MethylXGB requires exactly one tumor or tumor-mixture BAM. "
                  << "Use one -b input or add --disable-methyl-xgb for multi-BAM phasing.\n";
        die = true;
    }
    
    if ( opt::numThreads < 1 ){
        std::cerr << SUBPROGRAM " invalid threads. value: " 
                  << opt::numThreads 
                  << "\n please check -t, --threads=Num\n";
        die = true;
    }

    if ( opt::distance < 0 ){
        std::cerr << SUBPROGRAM " invalid distance. value: " 
                  << opt::distance 
                  << "\n please check -d or --distance=Num\n";
        die = true;
    }

    if ( opt::connectAdjacent < 0 ){
        std::cerr << SUBPROGRAM " invalid connectAdjacent. value: " 
                  << opt::connectAdjacent 
                  << "\n please check -a, --connectAdjacent=Num\n";
        die = true;
    }
    
    if ( opt::mappingQuality < 0 ){
        std::cerr << SUBPROGRAM " invalid mappingQuality. value: " 
                  << opt::mappingQuality 
                  << "\n please check -m, --mappingQuality=Num\n";
        die = true;
    }
    
    if ( opt::mismatchRate < 0 ){
        std::cerr << SUBPROGRAM " invalid mismatchRate. value: " 
                  << opt::mismatchRate 
                  << "\n please check -x, --mismatchRate=Num\n";
        die = true;
    }
    
    if ( opt::baseQuality < 0 ){
        std::cerr << SUBPROGRAM " invalid baseQuality. value: "
                  << opt::baseQuality
                  << "\n please check -m, --mappingQuality=[0~90]\n";
        die = true;
    }

    if ( opt::edgeWeight < 0 ){
        std::cerr << SUBPROGRAM " invalid edgeWeight. value: "
                  << opt::edgeWeight
                  << "\n please check -e, --edgeWeight=[0~1]\n";
        die = true;
    }

    if ( opt::baseQuality < 0 ){
        std::cerr << SUBPROGRAM " invalid baseQuality. value: "
                  << opt::baseQuality
                  << "\n please check -m, --mappingQuality=[0~90]\n";
        die = true;
    }

    if ( opt::edgeWeight < 0 ){
        std::cerr << SUBPROGRAM " invalid edgeWeight. value: "
                  << opt::edgeWeight
                  << "\n please check -m, --edgeWeight=[0~1]\n";
        die = true;
    }

    if ( opt::edgeThreshold < 0 || opt::edgeThreshold > 1 ){
        std::cerr << SUBPROGRAM " invalid edgeThreshold. value: " 
                  << opt::edgeThreshold 
                  << "\n please check -1, --edgeThreshold=[0~1]\n";
        die = true;
    }
    
    if ( opt::overlapThreshold < 0 || opt::overlapThreshold > 1 ){
        std::cerr << SUBPROGRAM " invalid overlapThreshold. value: " 
                  << opt::overlapThreshold 
                  << "\n please check -L, --overlapThreshold=[0~1]\n";
        die = true;
    }

    if ( opt::methylXgbExplicitEnable && opt::methylXgbExplicitDisable ){
        std::cerr << SUBPROGRAM ": --methyl-xgb and --disable-methyl-xgb cannot be used together.\n";
        die = true;
    }

    if ( opt::enableMethylXgb ){
        if( !std::isfinite(opt::methylXgbSnvThreshold) || opt::methylXgbSnvThreshold < 0.0 || opt::methylXgbSnvThreshold > 1.0 ){
            std::cerr << SUBPROGRAM " invalid methyl-xgb-snv-threshold. value: "
                      << opt::methylXgbSnvThreshold
                      << "\n please check --methyl-xgb-snv-threshold=[0~1]\n";
            die = true;
        }
        if( !std::isfinite(opt::methylXgbIndelThreshold) || opt::methylXgbIndelThreshold < 0.0 || opt::methylXgbIndelThreshold > 1.0 ){
            std::cerr << SUBPROGRAM " invalid methyl-xgb-indel-threshold. value: "
                      << opt::methylXgbIndelThreshold
                      << "\n please check --methyl-xgb-indel-threshold=[0~1]\n";
            die = true;
        }
        if( opt::methylXgbWindow <= 0 ){
            std::cerr << SUBPROGRAM " invalid methyl-window. value: "
                      << opt::methylXgbWindow
                      << "\n please check --methyl-window=Num\n";
            die = true;
        }
        if( !std::isfinite(opt::methylXgbMethLow) || !std::isfinite(opt::methylXgbMethHigh) || opt::methylXgbMethLow < 0.0 || opt::methylXgbMethHigh < 0.0 || opt::methylXgbMethLow >= opt::methylXgbMethHigh || opt::methylXgbMethHigh > 1.0 ){
            std::cerr << SUBPROGRAM " invalid methyl thresholds. meth-low: "
                      << opt::methylXgbMethLow
                      << " meth-high: "
                      << opt::methylXgbMethHigh
                      << "\n please check --meth-low and --meth-high\n";
            die = true;
        }
    }
    
    if ( opt::readConfidence < 0.5 || opt::readConfidence > 1 ){
        std::cerr << SUBPROGRAM " invalid readConfidence. value: " 
                  << opt::readConfidence 
                  << "\n please check -m, --readConfidence=[0.5~1]\n";
        die = true;
    }
    
    if ( opt::snpConfidence < 0.5 || opt::snpConfidence > 1 ){
        std::cerr << SUBPROGRAM " invalid snpConfidence. value: " 
                  << opt::snpConfidence 
                  << "\n please check -n, --snpConfidence=[0.5~1]\n";
        die = true;
    }
    
    if (opt::callerStr == "clairs_to_ss") {
        opt::caller = Caller::CLAIRS_TO_SS;
    } else if (opt::callerStr == "clairs_to_ssrs") {
        opt::caller = Caller::CLAIRS_TO_SSRS;
    } else if (opt::callerStr == "deepsomatic_to") {
        opt::caller = Caller::DEEPSOMATIC_TO;
    } else {
        std::cerr << SUBPROGRAM ": invalid caller option. Must be one of: clairs_to_ss, clairs_to_ssrs, deepsomatic_to\n";
        die = true;
    }

    if(opt::disableCalling){
        opt::somaticConnectAdjacent = 0;
    }

    if (die)
    {
        std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}

int PhasingMain(int argc, char** argv, std::string in_version)
{
    PhasingParameters ecParams;
    // set parameters
    PhasingOptions(argc, argv);
    // no file in command line

    ecParams.numThreads=opt::numThreads;
    ecParams.distance=opt::distance;
    ecParams.snpFile=opt::snpFile;
    ecParams.svFile=opt::svFile;
    ecParams.modFile=opt::modFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.generateDot=opt::generateDot;
    ecParams.phaseIndel=opt::phaseIndel;
    ecParams.caller=opt::caller;
    ecParams.callerStr=opt::callerStr;
    ecParams.disablePonTag=opt::disablePonTag;
    ecParams.disableCalling=opt::disableCalling;
    ecParams.disableRefineSomatic=opt::disableRefineSomatic;
    ecParams.enableMethylXgb=opt::enableMethylXgb;
    
    ecParams.connectAdjacent=opt::connectAdjacent;
    ecParams.mappingQuality=opt::mappingQuality;
    ecParams.mismatchRate=opt::mismatchRate;
    
    ecParams.baseQuality=opt::baseQuality;
    ecParams.edgeWeight=opt::edgeWeight;

    ecParams.edgeThreshold=opt::edgeThreshold;
    ecParams.overlapThreshold=opt::overlapThreshold;
    
    ecParams.snpConfidence=opt::snpConfidence;
    ecParams.readConfidence=opt::readConfidence;
    
    ecParams.ponFile=opt::ponFile;
    ecParams.strictPonFile=opt::strictPonFile;
    
    ecParams.somaticConnectAdjacent=opt::somaticConnectAdjacent;
    ecParams.methylXgbWindow=opt::methylXgbWindow;
    ecParams.methylXgbMethHigh=opt::methylXgbMethHigh;
    ecParams.methylXgbMethLow=opt::methylXgbMethLow;
    ecParams.methylXgbSnvThreshold=opt::methylXgbSnvThreshold;
    ecParams.methylXgbIndelThreshold=opt::methylXgbIndelThreshold;
    
    ecParams.version=in_version;
    ecParams.command=opt::command;
    
    ecParams.outputLOH = opt::outputLOH;
    ecParams.outputSGE = opt::outputSGE;
    ecParams.outputLGE = opt::outputLGE;
    ecParams.outputGE = opt::outputGE;
    ecParams.purity = opt::purity;

    PhasingProcess processor(ecParams);

    return 0;
}
