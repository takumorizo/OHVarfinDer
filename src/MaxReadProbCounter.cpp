//
//  MaxReadProbCounter.cpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/07/06.
//  Copyright © 2016年 森山卓也. All rights reserved.
//

#include "faidx.h"
#include "MutationCaller.h"
#include "BamReader.h"
#include "Parameters.h"
#include "MutationCallResult.h"
#include "CandidateWindow.h"
#include "HaplotypeBuilder.h"
#include "MutationModel.hpp"
#include "ProfileHMMUtils.h"
#include "ReadUtils.h"
#include "ReadsProfile.h"
#include "MaxReadProbCounter.h"
#include <algorithm>

MaxReadProbCounter::MaxReadProbCounter(BamReader &_tumorBamReader,
                         BamReader &_normalBamReader,
                         Parameters _params,
                         const faidx_t *_fai) : params(_params), tumorBamReader(_tumorBamReader),
normalBamReader(_normalBamReader), hapBuilder(_fai), fai(_fai) {
    hapBuilder.windowSize = 2 * (params.maxInsertSize + params.maxReadLength);
    //hapIndex;
    hapIndex.clear();
    //HapIndex O(+)H(-)
    hapIndex[std::make_pair(0, 0)] = std::make_pair(0, 0);
    hapIndex[std::make_pair(0, 1)] = std::make_pair(0, 1);
    hapIndex[std::make_pair(1, 0)] = std::make_pair(1, 0);
    hapIndex[std::make_pair(1, 1)] = std::make_pair(1, 1);
    hapIndex[std::make_pair(2, 0)] = std::make_pair(0, 0);
    hapIndex[std::make_pair(2, 1)] = std::make_pair(1, 1);
    hapIndex[std::make_pair(3, 0)] = std::make_pair(1, 0);
    hapIndex[std::make_pair(3, 1)] = std::make_pair(0, 1);

    //HapIndex O(-)H(+)
    hapIndex[std::make_pair(4, 0)] = std::make_pair(4, 0);
    hapIndex[std::make_pair(4, 1)] = std::make_pair(4, 1);
    hapIndex[std::make_pair(5, 0)] = std::make_pair(5, 0);
    hapIndex[std::make_pair(5, 1)] = std::make_pair(5, 1);
    hapIndex[std::make_pair(6, 0)] = std::make_pair(6, 0);
    hapIndex[std::make_pair(6, 1)] = std::make_pair(6, 1);
    hapIndex[std::make_pair(7, 0)] = std::make_pair(7, 0);
    hapIndex[std::make_pair(7, 1)] = std::make_pair(7, 1);

    //HapIndex O(+)H(+)
    hapIndex[std::make_pair(8, 0)] = std::make_pair(8, 0);
    hapIndex[std::make_pair(8, 1)] = std::make_pair(8, 1);
    hapIndex[std::make_pair(9, 0)] = std::make_pair(9, 0);
    hapIndex[std::make_pair(9, 1)] = std::make_pair(9, 1);
    hapIndex[std::make_pair(10, 0)] = std::make_pair(10, 0);
    hapIndex[std::make_pair(10, 1)] = std::make_pair(10, 1);
    hapIndex[std::make_pair(11, 0)] = std::make_pair(11, 0);
    hapIndex[std::make_pair(11, 1)] = std::make_pair(11, 1);
    hapIndex[std::make_pair(12, 0)] = std::make_pair(8, 0);
    hapIndex[std::make_pair(12, 1)] = std::make_pair(10, 1);
    hapIndex[std::make_pair(13, 0)] = std::make_pair(9, 0);
    hapIndex[std::make_pair(13, 1)] = std::make_pair(11, 1);
    hapIndex[std::make_pair(14, 0)] = std::make_pair(10, 0);
    hapIndex[std::make_pair(14, 1)] = std::make_pair(8, 1);
    hapIndex[std::make_pair(15, 0)] = std::make_pair(11, 0);
    hapIndex[std::make_pair(15, 1)] = std::make_pair(9, 1);
    
    //HapIndex O(-)H(-)S(-)
    hapIndex[std::make_pair(16, 0)] = std::make_pair(16, 0);
    hapIndex[std::make_pair(16, 1)] = std::make_pair(16, 1);
    hapIndex[std::make_pair(17, 0)] = std::make_pair(16, 0);
    hapIndex[std::make_pair(17, 1)] = std::make_pair(16, 1);
    hapIndex[std::make_pair(18, 0)] = std::make_pair(18, 0);
    hapIndex[std::make_pair(18, 1)] = std::make_pair(18, 1);

    //HapIndex O(-)H(-)S(-)
    hapIndex[std::make_pair(19, 0)] = std::make_pair(19, 0);
    hapIndex[std::make_pair(19, 1)] = std::make_pair(19, 1);
    hapIndex[std::make_pair(20, 0)] = std::make_pair(19, 0);
    hapIndex[std::make_pair(20, 1)] = std::make_pair(19, 1);
    hapIndex[std::make_pair(21, 0)] = std::make_pair(21, 0);
    hapIndex[std::make_pair(21, 1)] = std::make_pair(21, 1);
}

std::string MaxReadProbCounter::count(const CandidateWindow &window) {
    // prepare result
    MutationCallResult result(window);
    
    //
    // get target chr/pos from window
    //
    int pos = window.target.startInGenome;
    std::string chr = window.info.chr;
    LOG(logINFO) << chr << ":" << pos << " " << window.target << std::endl;
    LOG(logINFO) << "close variants: ";
    for (std::vector<Variant>::const_iterator it = window.closeVariants.begin(); it != window.closeVariants.end(); it++) {
        LOGP(logINFO) << *it << " ";
    }
    LOGP(logINFO) << std::endl;

    std::ostringstream tempOs;
    std::string ans;
    std::string posInfo;
    tempOs << chr << '\t' << window.target.startInGenome << '\t' << window.target.endInGenome << '\t' << window.target.ref << '\t' << window.target.obs << '\t';
    posInfo = tempOs.str();
    tempOs << "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
    tempOs << "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
    ans = tempOs.str();

    
    if (window.closeVariants.empty()) {
        result.setStatus("no_close_germline_variants");
    }
    bool ignoreBFComputation = (params.withoutBayesFactor || (window.info.ref_count_tumor == "-1") );
    
    std::vector<IRead *> allTumorReads, allNormalReads;
    
    try {
        //
        // get reads around the target
        // NOTE: these IRead * must be deleted manually.
        // adjast the window size by maxReads and variants
        // In OVerfinDer Calling, ignore SNPs nearby for fetching reads
        std::vector<Variant> closeVariantsForFetchingReads = window.closeVariants;
        
        if (ignoreBFComputation) {
            closeVariantsForFetchingReads.clear();
        }
        
        setTNReads(chr, pos, params.maxInsertSize, window, 
                   closeVariantsForFetchingReads, allTumorReads, allNormalReads, result);

        
        //
        // check reads
        //
        ReadsProfile allTumorReadsProfile(allTumorReads, window.target);
        ReadsProfile allNormalReadsProfile(allNormalReads, window.target);
        LOG(logINFO) << "profiles of all fetched reads:" << std::endl;
        LOGP(logINFO) << "tumor reads: " << allTumorReadsProfile << std::endl;
        LOGP(logINFO) << "normal reads: " << allNormalReadsProfile << std::endl;
        
        //
        // reads that cover the target position
        //
        std::vector<IRead *> tmpTumorReads = ReadUtils::coverPosition(allTumorReads, pos);
        std::vector<IRead *> tmpNormalReads = ReadUtils::coverPosition(allNormalReads, pos);
        
        //
        // check close variants before constructing haplotypes
        //
        bool isHaplotypeGoodState = checkHaplotypeGoodState(window, chr, tmpTumorReads, tmpNormalReads, params, result);
        //
        // construct local normal haplotypes
        //
        LOG(logINFO) << "constructing local normal haplotypes..." << std::endl;
        std::vector<Haplotype> normalHaps;
        std::vector<Haplotype> refHaps;
        if (not checkHaplotypeBuildPossibilityAndBuildHaplotype(normalHaps, refHaps, result, chr, closeVariantsForFetchingReads,
                                                                window, allTumorReads, allNormalReads) ){
            CandidateWindow new_window = window;
            new_window.closeVariants.clear();
            return count(new_window);
        }

        
        
        //
        // filter by mapping quality
        //
        std::vector<IRead *> tumorReads = ReadUtils::filterByMapQuality(tmpTumorReads, params.mapQualThreshold);
        std::vector<IRead *> normalReads = ReadUtils::filterByMapQuality(tmpNormalReads, params.mapQualThreshold);
        ReadsProfile filteredTumorReadsProfile(tumorReads, window.target);
        ReadsProfile filteredNormalReadsProfile(normalReads, window.target);
        LOG(logINFO) << "profiles of reads after mapping-quality filter:" << std::endl;
        LOGP(logINFO) << "tumor reads: " << filteredTumorReadsProfile << std::endl;
        LOGP(logINFO) << "normal reads: " << filteredNormalReadsProfile << std::endl;
        
        if (tumorReads.size() < params.minReads || normalReads.size() < params.minReads) {
            LOG(logERROR) << "too few number of reads: tumor = " << tumorReads.size()
            << ", normal = " << normalReads.size() << std::endl;
            LOG(logERROR) << "maybe, you should consider check mapQualThreshold" << std::endl;
            throw std::string("too_few_reads");
        }
        
        //
        // add a somatic haplotype and an error haplotype
        //
        LOG(logINFO) << "add somatic & error haplotypes..." << std::endl;
        std::vector<Haplotype> haps;
        std::vector<Haplotype> haps_ref;
        if(not checkHaplotypeBuildPossibilityAndBuildSomaticHaplotype(haps, haps_ref, normalHaps, refHaps,  
                                                                      window, closeVariantsForFetchingReads, chr,
                                                                      allTumorReads, allNormalReads)){
            CandidateWindow new_window = window;
            new_window.closeVariants.clear();
            return count(new_window);
        }

        
        std::vector<std::pair<Haplotype, Haplotype> > hapPairs;
        std::vector<std::pair<Haplotype, Haplotype> > hapPairsBasic;
        buildHapPaires(hapPairs, hapPairsBasic, haps, haps_ref);

        std::vector<std::vector<int> > tumorReadTypes  = MaxReadProbCounter::setReadTypes(closeVariantsForFetchingReads, window.target, tumorReads);
        std::vector<std::vector<int> > normalReadTypes = MaxReadProbCounter::setReadTypes(closeVariantsForFetchingReads ,window.target, normalReads);
        
        //
        // calculate alignment likelihoods of reads against haplotypes for OVarfinDer Model
        //
        LOG(logINFO) << "calculating alignment likelihoods..." << std::endl;
        std::vector<std::vector<double> > tumorLiks  = ProfileHMMUtils::calcAllLikelihoods(hapPairs, tumorReads,  tumorReadTypes, hapIndex);            
        std::vector<std::vector<double> > normalLiks = ProfileHMMUtils::calcAllLikelihoods(hapPairs, normalReads, normalReadTypes, hapIndex);

        std::vector<int > idx_h(4,0);   idx_h[0]=4;idx_h[1]=5;idx_h[2]=6;idx_h[3]=7;
        std::vector<int > idx_oh(4,0);  idx_oh[0]=8;idx_oh[1]=9;idx_oh[2]=10;idx_oh[3]=11;
        std::vector<int > idx_oh_all(8,0);  idx_oh_all[0]=8;idx_oh_all[1]=9;idx_oh_all[2]=10;idx_oh_all[3]=11;idx_oh_all[4]=12;idx_oh_all[5]=13;idx_oh_all[6]=14;idx_oh_all[7]=15;

        std::vector<int> typeBoth(22,0);
        typeBoth[8]  = typeBoth[9]  = typeBoth[10] = typeBoth[11] = typeBoth[12] = typeBoth[13] = typeBoth[14] = typeBoth[15] = 1;

        //
        // calculate Bayes factor, considering overlapping of paired-end read and HeteroSNP
        //
        try {
            HaplotypeBuilder::State state = HaplotypeBuilder::seeHaplotypeStateWithAlignmentCorrection(
                                                                                                 &normalLiks,
                                                                                                 &tumorLiks,
                                                                                                 normalReads,
                                                                                                 tumorReads,
                                                                                                 normalReadTypes,
                                                                                                 tumorReadTypes,
                                                                                                 window.target,
                                                                                                 &result.numHaplotypeInformativeReads,
                                                                                                 idx_h,idx_oh);
            LOG(logINFO) << "end check haplotype state " << std::endl;
            if( !isHaplotypeGoodState ){
                LOG(logERROR) << "something wrong state is already set" << std::endl;
            }else if( state == HaplotypeBuilder::SOMETHING_WRONG){
                LOG(logERROR) << "something wrong during calculation of bayes factor" << std::endl;
                throw std::string("something wrong during calculation of bayes factor");
            }else if(state == HaplotypeBuilder::LACK_HAPLOTYPE_INFO){
                result.setStatus("too_few_haplotype_informative_variant_supporting_reads");
                LOG(logERROR) << "too_few_haplotype_informative_variant_supporting_reads" << std::endl;
            }else if(state == HaplotypeBuilder::UMBALANCED_NORMAL_HAPLOTYPE){
                result.setStatus("umbalanced_normal_haplotypes");
                LOG(logERROR) << "umbalanced_normal_haplotypes" << std::endl;
            }else if( state == HaplotypeBuilder::KEEP_HAPLOTYPE or state == HaplotypeBuilder::SWITCH_HAPLOTYPE){
                if(state == HaplotypeBuilder::SWITCH_HAPLOTYPE) {
                   swapSNP(tumorLiks); swapSNP(normalLiks);
                   HaplotypeBuilder::swapHaps(hapPairs, idx_h, idx_oh_all);
                }
            }else{
                LOG(logINFO) << "Counting @ after seeHaplotypeStateWithAlignmentCorrection " << std::endl;

                // // Read Likelihood vote count start.
                // tumorReadTypes  = MaxReadProbCounter::setReadTypes(closeVariantsForFetchingReads, window.target, tumorReads,  true);
                // normalReadTypes = MaxReadProbCounter::setReadTypes(closeVariantsForFetchingReads, window.target, normalReads, true);
                // tumorLiks  = ProfileHMMUtils::calcAllLikelihoods(hapPairsBasic, tumorReads,  tumorReadTypes);
                // normalLiks = ProfileHMMUtils::calcAllLikelihoods(hapPairsBasic, normalReads, normalReadTypes);
            }

            // std::string getMaxReadProbSummary(std::vector<IRead *> reads, std::vector<std::vector<int> > readTypes);
            ans = posInfo;
            ans += getMaxReadProbSummary(tumorReads,  tumorReadTypes,  typeBoth, hapPairs);
            ans += getMaxReadProbSummary(normalReads, normalReadTypes, typeBoth, hapPairs);

            // {
            //     std::vector<IRead *>              tumorReadsOverlap;
            //     std::vector<std::vector<int> >    tumorReadTypesOverlap;
            //     std::vector<std::vector<double> > tumorlikesOverlap;

            //     std::vector<IRead *>              tumorReadsHap;
            //     std::vector<std::vector<int> >    tumorReadTypesHap;
            //     std::vector<std::vector<double> > tumorlikesHap;

            //     std::vector<IRead *>              tumorReadsBoth;
            //     std::vector<std::vector<int> >    tumorReadTypesBoth;
            //     std::vector<std::vector<double> > tumorlikesBoth;

            //     LOG(logINFO) << "Counting @ after init variables " << std::endl;

            //     collectOverlapReads(tumorReads, tumorReadTypes, tumorReadsOverlap);
            //     collectHapReads(tumorReads, tumorReadTypes, tumorReadsHap);
            //     collectOverlapHapReads(tumorReads, tumorReadTypes, tumorReadsBoth);
            //     LOG(logINFO) << "Counting @ after collection " << std::endl;
                
            //     LOG(logINFO) << " tumorReadsOverlap.size() " << tumorReadsOverlap.size() << ", tumorReadsHap.size() " << tumorReadsHap.size() << ", tumorReadsBoth.size() " << tumorReadsBoth.size() << std::endl;
                
            //     for(int i = 0; i < tumorReadsOverlap.size(); i++ ){ tumorReadTypesOverlap.push_back(typeBoth); }
            //     for(int i = 0; i < tumorReadsHap.size();     i++ ){ tumorReadTypesHap.push_back(typeBoth);     }
            //     for(int i = 0; i < tumorReadsBoth.size();    i++ ){ tumorReadTypesBoth.push_back(typeBoth);    }
            //     tumorlikesOverlap  = ProfileHMMUtils::calcAllLikelihoods(hapPairs, tumorReadsOverlap, tumorReadTypesOverlap, hapIndex);
            //     tumorlikesHap      = ProfileHMMUtils::calcAllLikelihoods(hapPairs, tumorReadsHap,     tumorReadTypesHap,     hapIndex);
            //     tumorlikesBoth     = ProfileHMMUtils::calcAllLikelihoods(hapPairs, tumorReadsBoth,    tumorReadTypesBoth,    hapIndex);

            //     LOG(logINFO) << " tumorlikesOverlap.size() " << tumorlikesOverlap.size() << ", tumorReadsHap.size() " << tumorlikesHap.size() << ", tumorlikesBoth.size() " << tumorlikesBoth.size() << std::endl;
            //     LOG(logINFO) << "Counting @ after alignment likelihoods computation " << std::endl;
            //     ans = posInfo;
            //     ans += maxReadProbCount(tumorlikesOverlap, tumorReadTypesOverlap, typeBoth);
            //     ans += maxReadProbCount(tumorlikesHap,     tumorReadTypesHap,     typeBoth);
            //     ans += maxReadProbCount(tumorlikesBoth,    tumorReadTypesBoth,    typeBoth);
            // }

            // {
            //     std::vector<IRead *>              normalReadsOverlap;
            //     std::vector<std::vector<int> >    normalReadTypesOverlap;
            //     std::vector<std::vector<double> > normallikesOverlap;

            //     std::vector<IRead *>              normalReadsHap;
            //     std::vector<std::vector<int> >    normalReadTypesHap;
            //     std::vector<std::vector<double> > normallikesHap;

            //     std::vector<IRead *>              normalReadsBoth;
            //     std::vector<std::vector<int> >    normalReadTypesBoth;
            //     std::vector<std::vector<double> > normallikesBoth;

            //     LOG(logINFO) << "Counting @ after init variables " << std::endl;
                
            //     collectOverlapReads(normalReads, normalReadTypes, normalReadsOverlap);
            //     collectHapReads(normalReads, normalReadTypes, normalReadsHap);
            //     collectOverlapHapReads(normalReads, normalReadTypes, normalReadsBoth);
            //     LOG(logINFO) << "Counting @ after collection " << std::endl;
            //     LOG(logINFO) << " normalReadsOverlap.size() " << normalReadsOverlap.size() << ", normalReadsHap.size() " << normalReadsHap.size() << ", normalReadsBoth.size() " << normalReadsBoth.size() << std::endl;
            //     for(int i = 0; i < normalReadsOverlap.size(); i++ ){ normalReadTypesOverlap.push_back(typeBoth); }
            //     for(int i = 0; i < normalReadsHap.size();     i++ ){ normalReadTypesHap.push_back(typeBoth);     }
            //     for(int i = 0; i < normalReadsBoth.size();    i++ ){ normalReadTypesBoth.push_back(typeBoth);    }
            //     normallikesOverlap  = ProfileHMMUtils::calcAllLikelihoods(hapPairs, normalReadsOverlap, normalReadTypesOverlap, hapIndex);
            //     normallikesHap      = ProfileHMMUtils::calcAllLikelihoods(hapPairs, normalReadsHap,     normalReadTypesHap,     hapIndex);
            //     normallikesBoth     = ProfileHMMUtils::calcAllLikelihoods(hapPairs, normalReadsBoth,    normalReadTypesBoth,    hapIndex);
            //     LOG(logINFO) << " normallikesOverlap.size() " << normallikesOverlap.size() << ", normalReadsHap.size() " << normallikesHap.size() << ", normallikesBoth.size() " << normallikesBoth.size() << std::endl;
            //     LOG(logINFO) << "Counting @ after alignment likelihoods computation " << std::endl;
            //     ans += maxReadProbCount(normallikesOverlap, normalReadTypesOverlap, typeBoth);
            //     ans += maxReadProbCount(normallikesHap,     normalReadTypesHap,     typeBoth);
            //     ans += maxReadProbCount(normallikesBoth,    normalReadTypesBoth,    typeBoth);
            // }

            LOG(logINFO) << "Counted" << std::endl;
        }catch (const std::exception& ex) {
            LOG(logERROR) << ex.what() << std::endl;
            result.setStatus("something_wrong", true);
        }catch (...) {
            LOG(logERROR) << "something wrong during read Counting" << std::endl;
            result.setStatus("something_wrong", true);
        }

    }catch (const std::exception& ex) {
        LOG(logERROR) << ex.what() << std::endl;
    }catch (std::string &s) {
        if (s == "too_many_reads_in_window") {
            result.setStatus(s, true);
        } else if (s == "too_few_reads") {
            result.setStatus(s, true);
        }else{
            result.setStatus("something_wrong", true);
        }
        LOG(logERROR) << s << std::endl;
        LOG(logINFO) << "UnCounted @ computing BF section" << std::endl;
        goto FINALLY; // C++ does not have the 'finally' statement.
    }
    
FINALLY:
    //
    // fin.
    //
    for (std::vector<IRead *>::iterator it = allTumorReads.begin();it != allTumorReads.end();it++) {
        delete *it;
    }
    for (std::vector<IRead *>::iterator it = allNormalReads.begin();it != allNormalReads.end();it++) {
        delete *it;
    }
    return ans.substr(0, std::max(0, (int)ans.length()-1)); // delete last \t
}

std::string MaxReadProbCounter::maxReadProbCount(const std::vector<std::vector<double> > &readLikes, const std::vector<std::vector<int> > &types, const std::vector<int> &maskType){
    std::vector<double> MaxReadProbCounted(22,0.0);
    for(int i = 0; i < readLikes.size(); i++){
        int    maxIndex  = -1;
        double maxLnProb = -1 * 1e300; 

        for(int j = 0; j < 22 ; j++)if(types[i][j] != 0)if(readLikes[i][j] > maxLnProb){    
            maxIndex = j;
            maxLnProb = readLikes[i][j]; 
        }
        int maxNum = 0;
        for(int j = 0; j < 22 ; j++)if(types[i][j] != 0)if(readLikes[i][j] == maxLnProb){
            maxNum++;
        }
        for(int j = 0; j < 22 ; j++)if(types[i][j] != 0)if(readLikes[i][j] == maxLnProb){
            MaxReadProbCounted[j] += 1.0 / maxNum;
        }
    }
    std::ostringstream ansOs;
    for(int j = 0; j < 22 ; j++)if(maskType[j] != 0){ 
        ansOs << MaxReadProbCounted[j] << '\t'; 
    }
    return ansOs.str();
}

std::string MaxReadProbCounter::getMaxReadProbSummary(std::vector<IRead *> reads, std::vector<std::vector<int> > readTypes, std::vector<int> typeBoth, std::vector<std::pair<Haplotype, Haplotype> > hapPairs){
    std::string ans = "";
    std::vector<IRead *>              readsOverlap;
    std::vector<std::vector<int> >    readTypesOverlap;
    std::vector<std::vector<double> > likesOverlap;

    std::vector<IRead *>              readsHap;
    std::vector<std::vector<int> >    readTypesHap;
    std::vector<std::vector<double> > likesHap;

    std::vector<IRead *>              readsBoth;
    std::vector<std::vector<int> >    readTypesBoth;
    std::vector<std::vector<double> > likesBoth;

    LOG(logINFO) << "Counting @ after init variables " << std::endl;

    collectOverlapReads(reads, readTypes, readsOverlap);
    collectHapReads(reads, readTypes, readsHap);
    collectOverlapHapReads(reads, readTypes, readsBoth);
    LOG(logINFO) << "Counting @ after collection " << std::endl;
    
    LOG(logINFO) << " tumorReadsOverlap.size() " << readsOverlap.size() << ", tumorReadsHap.size() " << readsHap.size() << ", tumorReadsBoth.size() " << readsBoth.size() << std::endl;
    
    for(int i = 0; i < readsOverlap.size(); i++ ){ readTypesOverlap.push_back(typeBoth); }
    for(int i = 0; i < readsHap.size();     i++ ){ readTypesHap.push_back(typeBoth);     }
    for(int i = 0; i < readsBoth.size();    i++ ){ readTypesBoth.push_back(typeBoth);    }
    likesOverlap  = ProfileHMMUtils::calcAllLikelihoods(hapPairs, readsOverlap, readTypesOverlap, hapIndex);
    likesHap      = ProfileHMMUtils::calcAllLikelihoods(hapPairs, readsHap,     readTypesHap,     hapIndex);
    likesBoth     = ProfileHMMUtils::calcAllLikelihoods(hapPairs, readsBoth,    readTypesBoth,    hapIndex);

    LOG(logINFO) << " tumorlikesOverlap.size() " << likesOverlap.size() << ", tumorReadsHap.size() " << likesHap.size() << ", tumorlikesBoth.size() " << likesBoth.size() << std::endl;
    LOG(logINFO) << "Counting @ after alignment likelihoods computation " << std::endl;
    ans += maxReadProbCount(likesOverlap, readTypesOverlap, typeBoth);
    ans += maxReadProbCount(likesHap,     readTypesHap,     typeBoth);
    ans += maxReadProbCount(likesBoth,    readTypesBoth,    typeBoth);
    return ans;
}

void MaxReadProbCounter::setTNReads(std::string chr, int pos,
                        int maxInsertSize, const CandidateWindow &window, std::vector<Variant> &closeVariants,
                        std::vector<IRead *> &allTumorReads, std::vector<IRead *> &allNormalReads,
                        MutationCallResult& result){
    while (true) {
        try {
            std::pair<int, int> windowRange = VariantUtils::getRange(window.target,
                                                                     closeVariants);
            LOG(logINFO) << "fetching tumor paired-reads..." << std::endl;            
            allTumorReads = tumorBamReader.getPairedReads(chr, pos,
                                                          windowRange.first,
                                                          windowRange.second,
                                                          params.maxInsertSize);
            LOG(logINFO) << "fetching normal paired-reads..." << std::endl;
            allNormalReads = normalBamReader.getPairedReads(chr, pos,
                                                            windowRange.first,
                                                            windowRange.second,
                                                            params.maxInsertSize);
            break;
        } catch (std::string &s) {
            if (s == "too_many_reads_in_window") {
                if (closeVariants.empty()) {
                    throw s;
                } else {
                    LOG(logWARNING) << "excluding furthest germline variants" << std::endl;
                    std::vector<Variant> tmp = VariantUtils::excludeFurthestVariant(window.target, closeVariants);
                    std::swap(tmp, closeVariants);
                    LOG(logWARNING) << "try fetching again" << std::endl;
                    continue;
                }
            } else {
                LOG(logERROR) << "error during fethcing reads: " << s << std::endl;
                // result.setBayesFactorBasic(-10000);
                throw s;
            }
        }
    }
}

bool MaxReadProbCounter::checkHaplotypeGoodState(const CandidateWindow &window, std::string chr, 
                             std::vector<IRead *>& tmpTumorReads, std::vector<IRead *>& tmpNormalReads, 
                             Parameters params, MutationCallResult &result){
    try {
        HaplotypeBuilder::checkCloseVariants(window.target,
                                             window.closeVariants,
                                             chr, tmpTumorReads,
                                             tmpNormalReads, params, result);
    } catch (std::string &s) {
        if (s == "no_close_germline_variants") {
            result.setStatus(s);
        } else if (s == "low_mapping_quality") {
            result.setStatus(s);
        } else if (s == "germline_variant_too_close") {
            result.setStatus(s);
        } else if (s == "germline_variants_overlapped") {
            result.setStatus(s);
        } else if (s == "germline_indel_too_close") {
            result.setStatus(s);
        } else if (s == "too_many_softclips_nearby") {
            result.setStatus(s);
        } else if (s == "too_many_indels_nearby") {
            result.setStatus(s);
        }
        return false;
    }
    return true;
}

bool MaxReadProbCounter::checkHaplotypeBuildPossibilityAndBuildHaplotype(std::vector<Haplotype> &normalHaps, std::vector<Haplotype> &refHaps, MutationCallResult &result,
                                                                   std::string chr, 
                                                                   const std::vector<Variant>& closeVariantsForFetchingReads, const CandidateWindow &window, 
                                                                   const std::vector<IRead *> &allTumorReads, const std::vector<IRead *> &allNormalReads){
    LOG(logINFO) << "constructing local normal haplotypes..." << std::endl;
    try {
        normalHaps = hapBuilder.constructNormalHaplotypes(window.target, closeVariantsForFetchingReads, chr,
                                                          allTumorReads, allNormalReads,
                                                          &result.MECScore, &result.MECnumReads);
        refHaps = hapBuilder.constructNormalHaplotypes(window.target, std::vector<Variant>(), chr,
                                                       allTumorReads, allNormalReads,
                                                       &result.MECScore, &result.MECnumReads);        
    } catch (std::string &s) {
        LOG(logERROR) << "something strange had happpend during custructing normal haplotypes" << std::endl;
        if (s == "error_making_haplotype") {
            LOG(logERROR) << "try again withoug close germline variants." << std::endl;
            return false;
        } else {
            throw s;
        }
    }
    return true;
}

bool MaxReadProbCounter::checkHaplotypeBuildPossibilityAndBuildSomaticHaplotype(std::vector<Haplotype> &haps, std::vector<Haplotype> &haps_ref,
                                                       const std::vector<Haplotype> &normalHaps,  const std::vector<Haplotype> &refHaps,  
                                                       const CandidateWindow &window, const std::vector<Variant> &closeVariantsForFetchingReads, std::string chr,
                                                       const std::vector<IRead *> &allTumorReads, const std::vector<IRead *> &allNormalReads){
    LOG(logINFO) << "add somatic & error haplotypes..." << std::endl;
    try {
        LOG(logINFO) << window.target << " : this is a target variant for adding haplotype " << std::endl;
        haps = hapBuilder.addSomaticAndErrorHaplotypes(normalHaps,
                                                       window.target, closeVariantsForFetchingReads, chr,
                                                       allTumorReads, allNormalReads);
        haps_ref = hapBuilder.addSomaticAndErrorHaplotypes(refHaps,
                                                           window.target, std::vector<Variant>(), chr,
                                                           allTumorReads, allNormalReads);
    } catch (std::string &s) {
        LOG(logERROR) << "something strange had happpend during custructing somatic and error haplotypes" << std::endl;
        if (s == "error_making_haplotype") {
            LOG(logERROR) << "try again withoug close germline variants." << std::endl;
            return false;
            // CandidateWindow new_window = window;
            // new_window.closeVariants.clear();
            // return call(new_window);
        } else {
            throw s;
        }
    }
    return true;
}

void MaxReadProbCounter::buildHapPaires(std::vector<std::pair<Haplotype, Haplotype> > &hapPairs, 
                                  std::vector<std::pair<Haplotype, Haplotype> > &hapPairsBasic,
                                  const std::vector<Haplotype> &hapsWithSNP, const std::vector<Haplotype> &hapsWithoutSNP){
    std::pair<Haplotype, Haplotype > hapRRefRef; hapRRefRef.first = hapsWithoutSNP[0]; hapRRefRef.second = hapsWithoutSNP[0];
    std::pair<Haplotype, Haplotype > hapRObsObs; hapRObsObs.first = hapsWithoutSNP[2]; hapRObsObs.second = hapsWithoutSNP[2];
    std::pair<Haplotype, Haplotype > hapRRefObs; hapRRefObs.first = hapsWithoutSNP[0]; hapRRefObs.second = hapsWithoutSNP[2];
    std::pair<Haplotype, Haplotype > hapRObsRef; hapRObsRef.first = hapsWithoutSNP[2]; hapRObsRef.second = hapsWithoutSNP[0];
    
    std::pair<Haplotype, Haplotype > hapARefRef; hapARefRef.first = hapsWithSNP[0]; hapARefRef.second = hapsWithSNP[0];
    std::pair<Haplotype, Haplotype > hapAObsObs; hapAObsObs.first = hapsWithSNP[2]; hapAObsObs.second = hapsWithSNP[2];
    std::pair<Haplotype, Haplotype > hapARefObs; hapARefObs.first = hapsWithSNP[0]; hapARefObs.second = hapsWithSNP[2];
    std::pair<Haplotype, Haplotype > hapAObsRef; hapAObsRef.first = hapsWithSNP[2]; hapAObsRef.second = hapsWithSNP[0];
    
    std::pair<Haplotype, Haplotype > hapBRefRef; hapBRefRef.first = hapsWithSNP[1]; hapBRefRef.second = hapsWithSNP[1];
    std::pair<Haplotype, Haplotype > hapBObsObs; hapBObsObs.first = hapsWithSNP[3]; hapBObsObs.second = hapsWithSNP[3];
    std::pair<Haplotype, Haplotype > hapBRefObs; hapBRefObs.first = hapsWithSNP[1]; hapBRefObs.second = hapsWithSNP[3];
    std::pair<Haplotype, Haplotype > hapBObsRef; hapBObsRef.first = hapsWithSNP[3]; hapBObsRef.second = hapsWithSNP[1];
    {
        // Overlap(+), HeteroSNP(-) pattern
        hapPairs.push_back(hapRRefRef);
        hapPairs.push_back(hapRObsObs);
        hapPairs.push_back(hapRRefObs);
        hapPairs.push_back(hapRObsRef);
        
        // Overlap(-), HeteroSNP(+) pattern
        hapPairs.push_back(hapARefRef);
        hapPairs.push_back(hapBRefRef);
        hapPairs.push_back(hapAObsObs);
        hapPairs.push_back(hapBObsObs);
        
        // Overlap(+), HeteroSNP(+) pattern
        hapPairs.push_back(hapARefRef); hapPairs.push_back(hapBRefRef);
        hapPairs.push_back(hapAObsObs); hapPairs.push_back(hapBObsObs);
        hapPairs.push_back(hapARefObs); hapPairs.push_back(hapBRefObs);
        hapPairs.push_back(hapAObsRef); hapPairs.push_back(hapBObsRef);
        
        // Overlap(-), HeteroSNP(-), plus obs, minus obs pattern
        hapPairs.push_back(hapRRefRef);
        hapPairs.push_back(hapRRefRef);
        hapPairs.push_back(hapRObsObs);

        hapPairs.push_back(hapRRefRef);
        hapPairs.push_back(hapRRefRef);
        hapPairs.push_back(hapRObsObs);        
    }
    {
        // Overlap(+), HeteroSNP(-) pattern
        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRObsObs);
        hapPairsBasic.push_back(hapRRefObs);
        hapPairsBasic.push_back(hapRObsRef);
        
        // Overlap(-), HeteroSNP(+) pattern
        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRObsObs);
        hapPairsBasic.push_back(hapRObsObs);
        
        // Overlap(+), HeteroSNP(+) pattern
        hapPairsBasic.push_back(hapARefRef); hapPairsBasic.push_back(hapBRefRef);
        hapPairsBasic.push_back(hapAObsObs); hapPairsBasic.push_back(hapBObsObs);
        hapPairsBasic.push_back(hapARefObs); hapPairsBasic.push_back(hapBRefObs);
        hapPairsBasic.push_back(hapAObsRef); hapPairsBasic.push_back(hapBObsRef);
        
        // Overlap(-), HeteroSNP(-), plus obs, minus obs pattern
        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRObsObs);

        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRRefRef);
        hapPairsBasic.push_back(hapRObsObs);        
    }
}


void MaxReadProbCounter::collectOverlapReads(const std::vector<IRead *> &reads, const std::vector<std::vector<int> > &types, std::vector<IRead *> &readsOverlap){
    readsOverlap.clear();
    for(int i = 0; i < reads.size(); i++)if(types[i][0] != 0){
        readsOverlap.push_back(reads[i]);
    }
}

void MaxReadProbCounter::collectHapReads(const std::vector<IRead *> &reads, const std::vector<std::vector<int> > &types, std::vector<IRead *> &readsHap){
    readsHap.clear();
    for(int i = 0; i < reads.size(); i++)if(types[i][4] != 0){
        readsHap.push_back(reads[i]);
    }
}

void MaxReadProbCounter::collectOverlapHapReads(const std::vector<IRead *> &reads, const std::vector<std::vector<int> > &types, std::vector<IRead *> &readsBoth){
    readsBoth.clear();
    for(int i = 0; i < reads.size(); i++)if(types[i][8] != 0){
        readsBoth.push_back(reads[i]);
    }
}

MutationCallResult MaxReadProbCounter::call(const CandidateWindow &window) {
    // prepare result
    MutationCallResult result(window);
    return result;
}