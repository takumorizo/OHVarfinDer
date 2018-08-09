#ifndef ReadTest_Parameters_h
#define ReadTest_Parameters_h

#include <cstdlib>
#include <string>
#include <vector>
#include "log.h"
#include "Utils.h"
#include "cmdline.h"

class Parameters {
public:
    class BayesEMParameters { // for Bayesian EM algorithm
    public:
        std::vector<double> mut_a0, mut_b0, mut_c0, mut_d0, mut_e0, mut_f0, mut_g0, mut_h0, mut_i0, mut_j0;
        std::vector<double> err_a0, err_b0, err_c0, err_d0, err_e0, err_f0, err_g0, err_h0, err_i0, err_j0;
        BayesEMParameters()  {};
        ~BayesEMParameters() {};
        void clear_all() {
            mut_a0.clear(); mut_b0.clear();
            mut_c0.clear(); mut_d0.clear();
            mut_e0.clear(); mut_f0.clear();
            mut_g0.clear(); mut_h0.clear();
            mut_i0.clear(); mut_j0.clear();
            
            err_a0.clear(); err_b0.clear();
            err_c0.clear(); err_d0.clear();
            err_e0.clear(); err_f0.clear();
            err_g0.clear(); err_h0.clear();
            err_i0.clear(); err_j0.clear();
        }

        template <class X> void print_veci(std::vector<X> vec) const {
            for (int i = 0; i < vec.size(); i++) {
                LOGP(logINFO) << vec[i] << " ";
            }
        }

        void print() const{
            LOG(logINFO) << "\tMutationModel:" << std::endl;
            LOG(logINFO) << "\t\ta0: ";
            print_veci(mut_a0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tb0: ";
            print_veci(mut_b0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tc0: ";
            print_veci(mut_c0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\tErrorModel:" << std::endl;
            LOG(logINFO) << "\t\ta0: ";
            print_veci(err_a0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tb0: ";
            print_veci(err_b0);
            LOGP(logINFO) << std::endl;
            LOG(logINFO) << "\t\tc0: ";
            print_veci(err_c0);
            LOGP(logINFO) << std::endl;
        }

        BayesEMParameters &operator=(const BayesEMParameters &t) {
            mut_a0 = t.mut_a0;  mut_b0 = t.mut_b0;
            mut_c0 = t.mut_c0;  mut_d0 = t.mut_d0;
            mut_e0 = t.mut_e0;  mut_f0 = t.mut_f0;
            mut_g0 = t.mut_g0;  mut_h0 = t.mut_h0;
            mut_i0 = t.mut_i0;  mut_j0 = t.mut_j0;
            
            err_a0 = t.err_a0;  err_b0 = t.err_b0;
            err_c0 = t.err_c0;  err_d0 = t.err_d0;
            err_e0 = t.err_e0;  err_f0 = t.err_f0;
            err_g0 = t.err_g0;  err_h0 = t.err_h0;
            err_i0 = t.err_i0;  err_i0 = t.err_j0;

            return *this;
        }
    };
    class PileUpParameters{
    public:
        typedef enum {TUMOR,NORMAL,HETEROSNP,TRIALLELE} SAMPLE;
        SAMPLE sample;
        PileUpParameters(){
            init(TUMOR);
        }
        PileUpParameters(SAMPLE sample){
            init(sample);
        }
        void init(SAMPLE sample){
            this->sample = sample;
            if(sample == TUMOR){
                minBQ = 15;
                minDepth = 100; maxDepth = 1000000;
                minDepthPlus=0; maxDepthPlus=1000000; minDepthMinus=0; maxDepthMinus=1000000;
                minObsRate=0.01; maxObsRate=0.07;
                minObsRatePlus=0.0; maxObsRatePlus=1.0; minObsRateMinus=0.0; maxObsRateMinus=1.0;
                minObsNum=3; maxObsNum=1000000;
                minObsNumPlus=0; maxObsNumPlus=1000000; minObsNumMinus=0; maxObsNumMinus=1000000;
                minRefNum=0; maxRefNum=1000000;
                minRefNumPlus=0; maxRefNumPlus=1000000; minRefNumMinus=0; maxRefNumMinus=1000000;
                avgBaseQualityThreshold = 20;
            }else if(sample == NORMAL){
                minBQ = 15;
                minDepth = 100; maxDepth = 1000000;
                minDepthPlus=0; maxDepthPlus=1000000; minDepthMinus=0; maxDepthMinus=1000000;
                minObsRate=0.0; maxObsRate=0.01;
                minObsRatePlus=0.0; maxObsRatePlus=1.0; minObsRateMinus=0.0; maxObsRateMinus=1.0;
                minObsNum=0; maxObsNum=1;
                minObsNumPlus=0; maxObsNumPlus=1000000; minObsNumMinus=0; maxObsNumMinus=1000000;
                minRefNum=0; maxRefNum=1000000;
                minRefNumPlus=0; maxRefNumPlus=1000000; minRefNumMinus=0; maxRefNumMinus=1000000;
                avgBaseQualityThreshold = 20;
            }else if(sample == HETEROSNP){
                minBQ = 15;
                minDepth = 20; maxDepth = 1000000;
                minDepthPlus=0; maxDepthPlus=1000000; minDepthMinus=0; maxDepthMinus=1000000;
                minObsRate=0.2; maxObsRate=1.0;
                minObsRatePlus=0.0; maxObsRatePlus=1.0; minObsRateMinus=0.0; maxObsRateMinus=1.0;
                minObsNum=5; maxObsNum=1000000;
                minObsNumPlus=0; maxObsNumPlus=1000000; minObsNumMinus=0; maxObsNumMinus=1000000;
                minRefNum=0; maxRefNum=1000000;
                minRefNumPlus=0; maxRefNumPlus=1000000; minRefNumMinus=0; maxRefNumMinus=1000000;
                avgBaseQualityThreshold = 20;
            }else if(sample == TRIALLELE){
                minBQ = 15;
                minDepth = 20; maxDepth = 1000000;
                minDepthPlus=0; maxDepthPlus=1000000; minDepthMinus=0; maxDepthMinus=1000000;
                minObsRate=0.03; maxObsRate=1.0;
                minObsRatePlus=0.0; maxObsRatePlus=1.0; minObsRateMinus=0.0; maxObsRateMinus=1.0;
                minObsNum=4; maxObsNum=1000000;
                minObsNumPlus=0; maxObsNumPlus=1000000; minObsNumMinus=0; maxObsNumMinus=1000000;
                minRefNum=0; maxRefNum=1000000;
                minRefNumPlus=0; maxRefNumPlus=1000000; minRefNumMinus=0; maxRefNumMinus=1000000;
                avgBaseQualityThreshold = 20;
            }
        }
        int minBQ;
        int minDepth,maxDepth;
        int minDepthPlus, maxDepthPlus,minDepthMinus, maxDepthMinus;
        double minObsRate,maxObsRate;
        double minObsRatePlus,maxObsRatePlus,minObsRateMinus,maxObsRateMinus;
        int minObsNum, maxObsNum;
        int minObsNumPlus,maxObsNumPlus,minObsNumMinus,maxObsNumMinus;
        int minRefNum, maxRefNum;
        int minRefNumPlus,maxRefNumPlus,minRefNumMinus,maxRefNumMinus;
        double avgBaseQualityThreshold;
    };
    PileUpParameters pileN, pileT, pileHetero, pileTriAllele;
    
    // typedef enum { OVAR,OVARFINDER,OHVARFINDER,OHVARFINDER,HAPMUC,HETEROSNPCALL } ALGO;
    typedef enum { OHVARFINDER,HAPMUC,HETEROSNPCALL} ALGO;
    ALGO method;

    Parameters(int argc, const char *argv[]);
    
    void getFromCommandLineArguments(int argc, const char *argv[]);
    static void parseHyperParameters(std::vector<double> &vec, std::string s);
    
    // === shared parameters between at least two methods ===
    int maxReads, maxReadLength, minReads;
    int maxInsertSize, maxIndelLength;
    int minDistanceSNP, minDistanceGermlineIndel;
    double priorDelByError;
    double averageMapQualThreshold, mapQualThreshold, priorIndel, priorSNP, EMtol;
    int averageMapPhredQualThreshold, mapPhredQualThreshold;

    double pError, pMut;
    std::string outFilePrefix, refFileName, tumorBam, normalBam, windowFile;
    bool showHapAlignments, quiet;
    bool singleReads;
    double indelPosessionFreqThreshold, softClipPosessionFreqThreshold;
    bool withoutBayesFactor;

    std::string region,pileupFile;
    double heteroSNPConfidenceInterval;
    int pileUpBufferSize;


    int fReadFilter;
    int FReadFilter;
    
    // === parameters only for OHVarfinDer method ===
    BayesEMParameters OHVarfinDParams2;
    
    // === parameters only for HapMuc method ===
    BayesEMParameters hap3_params, hap2_params;
    
    // === parameters only for OHVarfinDParams2NoCategory method ===
    BayesEMParameters OHVarfinDParams2NoCategory;
private:
    void getFromCommandLineArgumentsInOHVarfinDer(           cmdline::parser& a,int argc, const char *argv[]);
    void getFromCommandLineArgumentsInHapMuC(                 cmdline::parser& a,int argc, const char *argv[]);
    void getFromCommandLineArgumentsInOHVarfinDerNoCategory( cmdline::parser& a,int argc, const char *argv[]);
};
#endif
