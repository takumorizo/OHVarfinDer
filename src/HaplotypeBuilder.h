/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#ifndef __HaplotypeBuilder__
#define __HaplotypeBuilder__

#include <vector>
#include <string>

#include "faidx.h"

#include "Haplotype.hpp"
#include "Variant.h"
#include "BasicHaplotypeBuilder.h"
#include "IRead.h"
#include "Parameters.h"
#include "MutationCallResult.h"

class HaplotypeBuilder : public BasicHaplotypeBuilder {
    const faidx_t *fai;
public:
    typedef enum { SOMETHING_WRONG, LACK_HAPLOTYPE_INFO, UMBALANCED_NORMAL_HAPLOTYPE, SWITCH_HAPLOTYPE,KEEP_HAPLOTYPE } State;

    int windowSize;
    std::vector<Haplotype> buildHaplotypes(const Variant &targetVariant, const std::vector<Variant> &closeVariants,
                                           std::string chr,
                                           const std::vector<IRead *> &tumorReads,
                                           const std::vector<IRead *> &normalReads);
    std::vector<Haplotype> constructNormalHaplotypes(const Variant &targetVariant,
                                                     const std::vector<Variant> &closeVariants,
                                                     std::string chr,
                                                     const std::vector<IRead *> &tumorReads,
                                                     const std::vector<IRead *> &normalReads,
                                                     int *_MECScore, int *_MECNumReads);
    std::vector<Haplotype> addSomaticAndErrorHaplotypes(const std::vector<Haplotype> &normalHaps,
                                                        const Variant &targetVariant,
                                                        const std::vector<Variant> &closeVariants,
                                                        std::string chr,
                                                        const std::vector<IRead *> &tumorReads,
                                                        const std::vector<IRead *> &normalReads);
    static std::vector<Variant> checkCloseVariants(const Variant &targetVariant,
                                            const std::vector<Variant> &closeVariants,
                                            std::string chr,
                                            const std::vector<IRead *> &tumorReads,
                                            const std::vector<IRead *> &normalReads,
                                            Parameters &params, MutationCallResult &result);
    HaplotypeBuilder(const faidx_t *fai);
    static void checkAndSort(std::vector<Haplotype> *haps,
                             std::vector<std::vector<double> > *normalLiks,
                             std::vector<std::vector<double> > *tumorLiks,
                             const std::vector<IRead *> &normalReads,
                             const std::vector<IRead *> &tumorReads,
                             const Variant &targetVariant,
                             int *numHaplotypeInfomativeReads);
    static bool checkSwap(std::vector<Haplotype> *haps,
                                        std::vector<std::vector<double> > *normalLiks,
                                        std::vector<std::vector<double> > *tumorLiks,
                                        const std::vector<IRead *> &normalReads,
                                        const std::vector<IRead *> &tumorReads,
                                        const Variant &targetVariant,
                                        int *numHaplotypeInfomativeReads,
                                        std::vector<int> idx);

    static State seeHaplotypeStateWithAlignmentCorrection(std::vector<std::vector<double> > *normalLiks,
                                                          std::vector<std::vector<double> > *tumorLiks,
                                                          const std::vector<IRead *> &normalReads,
                                                          const std::vector<IRead *> &tumorReads,
                                                          const std::vector<std::vector<int> > &normalReadTypes,
                                                          const std::vector<std::vector<int> > &tumorReadTypes,
                                                          const Variant &targetVariant,
                                                          int *numHaplotypeInfomativeReads,
                                                          const std::vector<int> &idx_h, const std::vector<int> &idx_oh);
    static std::vector<std::vector<double> > convertLiks2basic(const std::vector<std::vector<double> > &liks);
    static std::vector<Haplotype> convertHaps2basic(const std::vector<Haplotype> &haps);
    static void swapHaps(std::vector<std::pair<Haplotype, Haplotype> > &haps, std::vector<int> idx_h, std::vector<int> idx_oh);
    std::vector<Haplotype> prepareHaplotypesForBasic(std::string chr,const Variant &targetVariant);
private:
    static void correctLikelihood(std::vector<std::vector<double> > *liks, std::vector<std::vector<int> > types, double maxLikelihood, const std::vector<int> idx ){
        for(int i = 0; i < liks->size(); i++)if(types[i][idx[0]] != 0){
            bool flag = false;
            double max = liks->at(i)[idx[0]];
            for(int j = 1; j < 4; j++) max = std::max(max, liks->at(i)[idx[j]]);
            if( (max - maxLikelihood ) < -25.0 ) flag = true;
            if(flag){
                LOG(logDEBUG) << "filter low-quality reads" << std::endl;
                for (int j = 0;j < 4;j++) {
                    LOGP(logDEBUG) << liks->at(i)[idx[j]] << " ";
                    liks->at(i)[idx[j]] = max;
                }
                LOGP(logDEBUG) << std::endl;
            }
        }
    }
    static int checkHaplotypeInformativeReadsNum(std::vector<std::vector<double> > *liks, std::vector<std::vector<int> > types, const std::vector<int> &idx ){
        int count = 0;
        for(int i = 0; i < liks->size(); i++) if(types[i][idx[0]] != 0){
            double max = liks->at(i)[idx[0]];
            int max_index = 0;
            for(int j = 1; j < 4; j++){ if(max < liks->at(i)[idx[j]]){
                max = std::max(max, liks->at(i)[idx[j]]);
                max_index = j;
            }}
            if(std::abs(liks->at(i)[idx[2]]-liks->at(i)[idx[3]]) > 10e-2 && (max_index == 2 || max_index == 3) ){
                count++;
            }
        }
        return count;
    }

    static void flattenLowQualityReadLikelihoods(std::vector<std::vector<double> > *liks, std::vector<std::vector<int> > types, const std::vector<int> &idx, double maxLik, double diffMax){
        diffMax = std::abs(diffMax);
        for(int i = 0; i < liks->size(); i++) if(types[i][idx[0]] != 0){
            double max = liks->at(i)[idx[0]];
            int max_index = 0;
            for(int j = 1; j < 4; j++) if(max < liks->at(i)[idx[j]]){
                max = std::max(max, liks->at(i)[idx[j]]); 
                max_index = j;
            }
            if( maxLik - max > diffMax)for(int j = 0; j < 4; j++){
              liks->at(i)[idx[0]] = max;
            }
        }
    }

};

#endif /* defined(__HaplotypeBuilder__) */
