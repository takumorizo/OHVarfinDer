//
//  OHVarfinDer2.hpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/07/06.
//  Copyright © 2016年 森山卓也. All rights reserved.
//

#ifndef OHVarfinDer2_hpp
#define OHVarfinDer2_hpp

#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include "PairedRead.h"
#include "SingleRead.h"
#include "SamRead.h"
#include "ReadUtils.h"
#include "BaseMutationCaller.h"
#include "CandidateWindow.h"
#include "Variant.h"
#include "log.h"
#include <iostream>
#include <sstream>
#include <string>


class OHVarfinDer2 : public BaseMutationCaller{
public:
    HaplotypeBuilder hapBuilder;
    BamReader &tumorBamReader;
    BamReader &normalBamReader;
    const Parameters params;
    const faidx_t *fai;
    OHVarfinDer2(BamReader &tumorBamReader, BamReader &normalBamReader, Parameters params, const faidx_t *fai);
    virtual MutationCallResult call(const CandidateWindow &candidateWindow);

    // TODO : rewrite to the normal private function
    static inline std::vector<std::vector<int> > setReadTypes(const std::vector<Variant> snps,
                                                              const Variant &variant,
                                                              const std::vector<IRead*> &reads,
                                                              bool rmHeteroSNP = false) {
        
        // init ans
        std::vector<std::vector<int> > ans(reads.size());
        std::vector<int> type1(22,0); // overlap(+) heterosnp(-)
        std::vector<int> type2(22,0); // overlap(-) heterosnp(+)
        std::vector<int> type3(22,0); // overlap(+) heterosnp(+)
        std::vector<int> type4(22,0); // overlap(-) heterosnp(-), on Plus  Strand
        std::vector<int> type5(22,0); // overlap(-) heterosnp(-), on Minus Strand
        std::vector<int> type6(22,0); // non of above.
        
        type1[0]  = type1[1]  = type1[2]  = type1[3]  = 1;
        type2[4]  = type2[5]  = type2[6]  = type2[7]  = 1;
        type3[8]  = type3[9]  = type3[10] = type3[11] = type3[12] = type3[13] = type3[14] = type3[15] = 1;
        type4[16] = type4[17] = type4[18] = 1;
        type5[19] = type5[20] = type5[21] = 1;

        int t1=0,t2=0,t3=0,t4=0,t5=0,t6=0;
        
        for(int i = 0; i < reads.size(); i++){
            IRead* read = reads[i];
            bool onSNP = ReadUtils::coverVariants(read ,snps);
            bool overlap = false;
            bool onPlus  = false;
            bool onMinus = false;

            if(typeid(*read) == typeid(PairedRead) ){
                const PairedRead* pairRead = dynamic_cast<const PairedRead *>(read);
                SamRead* read1 = pairRead->first;   // plus  strand read
                SamRead* read2 = pairRead->second;  // minus strand read
                
                if( read1 != NULL && !read1->onReverseStrand) onPlus  = read1->cover(variant.startInGenome);
                if( read1 != NULL &&  read1->onReverseStrand) onMinus = read1->cover(variant.startInGenome);

                if( read2 != NULL && !read2->onReverseStrand) onPlus  = read2->cover(variant.startInGenome);
                if( read2 != NULL &&  read2->onReverseStrand) onMinus = read2->cover(variant.startInGenome);
                overlap = onPlus && onMinus;

            }else if(typeid(*read) == typeid(SingleRead)){
                 const SingleRead* singleRead = dynamic_cast<const SingleRead *>(read);
                SamRead* readSam = singleRead->read;
                if( readSam != NULL && !readSam->onReverseStrand) onPlus  = readSam->cover(variant.startInGenome);
                if( readSam != NULL &&  readSam->onReverseStrand) onMinus = readSam->cover(variant.startInGenome);

            }else{
                throw std::string("setReadTypes @ OVarfindDer : unexpected class type of IRead");
            }

            // if(         overlap && !onSNP ){ t2++; ans[i] = type2;}
            // else if(   !overlap &&  onSNP ){ t2++; ans[i] = type2;}
            // else if(    overlap &&  onSNP ){ t2++; ans[i] = type2;}
            // else if(   !overlap && !onSNP ){ t2++; ans[i] = type2;}

/*
            if(         overlap && !onSNP ){ t1++; ans[i] = type1;}
            else if(   !overlap &&  onSNP ){ t2++; ans[i] = type2;}
            else if(    overlap &&  onSNP ){ t3++; ans[i] = type3;}
            else if(   !overlap && !onSNP ){ t2++; ans[i] = type2;}
*/

            if (     !rmHeteroSNP ){
                if(         overlap && !onSNP ){ t1++; ans[i] = type1;}
                else if(   !overlap &&  onSNP ){ t2++; ans[i] = type2;}
                else if(    overlap &&  onSNP ){ t3++; ans[i] = type3;}
                else if(   !overlap && !onSNP ){ 
                    if(      onPlus){ t4++; ans[i] = type4;}
                    else if(onMinus){ t5++; ans[i] = type5;}
                    else            { t6++; ans[i] = type6;}
                }
            }else if( rmHeteroSNP ){
                if(        overlap && !onSNP ){ 
                    t1++; ans[i] = type1;
                }else if(  !overlap &&  onSNP ){ 
                    if(      onPlus){ t4++; ans[i] = type4;}
                    else if(onMinus){ t5++; ans[i] = type5;}
                    else            { t6++; ans[i] = type6;}
                }else if(   overlap &&  onSNP ){
                    t1++; ans[i] = type1;
                }else if(  !overlap && !onSNP ){
                    if(      onPlus){ t4++; ans[i] = type4;}
                    else if(onMinus){ t5++; ans[i] = type5;}
                    else            { t6++; ans[i] = type6;}
                }
            }
        }
        
        if(rmHeteroSNP){
            LOG(logINFO) << "ReadTypeProfile @ rmHeteroSNP " << std::endl;
            LOG(logINFO) << "Overlap - HeteroSNP - (+): " <<  t4 << std::endl;
            LOG(logINFO) << "Overlap - HeteroSNP - (-): " <<  t5 << std::endl;
            LOG(logINFO) << "Overlap + HeteroSNP -    : " <<  t1 << std::endl;
            LOG(logINFO) << "Overlap - HeteroSNP +    : " <<  t2 << std::endl;
            LOG(logINFO) << "Overlap + HeteroSNP +    : " <<  t3 << std::endl;
            LOG(logINFO) << "non of above             : " <<  t6 << std::endl;
        }else {
            LOG(logINFO) << " ReadTypeProfile @ not rmHeteroSNP " << std::endl;
            LOG(logINFO) << "Overlap - HeteroSNP - (+): " <<  t4 << std::endl;
            LOG(logINFO) << "Overlap - HeteroSNP - (-): " <<  t5 << std::endl;
            LOG(logINFO) << "Overlap + HeteroSNP -    : " <<  t1 << std::endl;
            LOG(logINFO) << "Overlap - HeteroSNP +    : " <<  t2 << std::endl;
            LOG(logINFO) << "Overlap + HeteroSNP +    : " <<  t3 << std::endl;
            LOG(logINFO) << "non of above             : " <<  t6 << std::endl;
        }
        
        return ans;
    }
private:
    const double fromLnToLog10Factor;
    std::map<std::pair<int, int>, std::pair<int, int> > hapIndex;
    void swapSNP(std::vector<std::vector<double> > &readLiks ){
         for(int i = 0; i < readLiks.size(); i++){
            // 4<->5, 6<->7,
            std::swap(readLiks[i][4],readLiks[i][5]);
            std::swap(readLiks[i][6],readLiks[i][7]);

            // 8<->9, 10<->11, 12<->13, 14<->15
            std::swap(readLiks[i][8],readLiks[i][9]);
            std::swap(readLiks[i][10],readLiks[i][11]);
            std::swap(readLiks[i][12],readLiks[i][13]);
            std::swap(readLiks[i][14],readLiks[i][15]);
         }
    }
    
    double getPriorRatio(CandidateWindow window){
        bool isIndel = ( window.target.ref == "-" || window.target.obs == "-" );
        if(isIndel){ return params.priorIndel;   }
        else {       return params.priorSNP;     }
    }

    void setTNReads(std::string chr, int pos,
                    int maxInsertSize, const CandidateWindow &window, std::vector<Variant> &closeVariants,
                    std::vector<IRead *> &allTumorReads, std::vector<IRead *> &allNormalReads,
                    MutationCallResult& result);

    bool checkHaplotypeGoodState(const CandidateWindow &window, std::string chr, 
                                 std::vector<IRead *>& tmpTumorReads, std::vector<IRead *>& tmpNormalReads, 
                                 Parameters params, MutationCallResult &result);

    bool checkHaplotypeBuildPossibilityAndBuildHaplotype(std::vector<Haplotype> &normalHaps, std::vector<Haplotype> &refHaps, MutationCallResult &result,
                                                         std::string chr, 
                                                         const std::vector<Variant>& closeVariantsForFetchingReads, const CandidateWindow &window, 
                                                         const std::vector<IRead *> &allTumorReads, const std::vector<IRead *> &allNormalReads);

    bool checkHaplotypeBuildPossibilityAndBuildSomaticHaplotype(std::vector<Haplotype> &haps, std::vector<Haplotype> &haps_ref,
                                                                const std::vector<Haplotype> &normalHaps,  const std::vector<Haplotype> &refHaps,  
                                                                const CandidateWindow &window, const std::vector<Variant> &closeVariantsForFetchingReads, std::string chr,
                                                                const std::vector<IRead *> &allTumorReads, const std::vector<IRead *> &allNormalReads);

    void buildHapPaires(std::vector<std::pair<Haplotype, Haplotype> > &hapPairs, 
                        std::vector<std::pair<Haplotype, Haplotype> > &hapPairsBasic,
                        const std::vector<Haplotype> &hapsWithSNP, const std::vector<Haplotype> &hapsWithoutSNP);

    bool isBFWithHapComputable(bool isHaplotypeGoodState, HaplotypeBuilder::State state, MutationCallResult &result);

    static double getFactorFromLnToLog10();
};
#endif /* OHVarfinDer2_hpp */
