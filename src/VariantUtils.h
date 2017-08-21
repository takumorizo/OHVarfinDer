#ifndef __haplotype_builder__VariantUtils__
#define __haplotype_builder__VariantUtils__

#include <vector>
#include <sstream>
#include <algorithm>

#include "Variant.h"
#include "VariantFromSAM.h"
#include "log.h"
#include "CandidateVariant.h"
#include <cstdlib>
#include "CandidateWindow.h"

class VariantUtils {
public:
    static inline Variant variantFromSAMToVariant(const VariantFromSAM &_v) {
        Variant v;
        v.startInGenome = _v.left;
        v.endInGenome = _v.right;
        v.length = _v.length;
        v.ref = _v.ref;
        v.obs = _v.obs;
        std::stringstream ss;
        if (_v.type == VariantFromSAM::MISMATCH) {
            ss << _v.ref << "=>" << _v.obs;
            v.type = Variant::SNP;
        } else if (_v.type == VariantFromSAM::INS) {
            ss << "+" << _v.obs;
            v.type = Variant::INS;
        } else if (_v.type == VariantFromSAM::DEL) {
            ss << "-" << _v.ref;
            v.type = Variant::DEL;
        } else {
            throw std::string("fromVariantFromSAM_not_supported");
        }
        v.seq = ss.str();
        v.originalString = v.seq;
        return v;
    }
    
    static inline VariantFromSAM variantToVariantFromSAM(const Variant &_v) {
        VariantFromSAM v;
        v.left = _v.startInGenome;
        v.right = _v.endInGenome;
        v.length = _v.length;
        v.ref = _v.ref;
        v.obs = _v.obs;
        if (_v.type == Variant::SNP) {
            v.type = VariantFromSAM::MISMATCH;
        } else if (_v.type == Variant::INS) {
            v.type = VariantFromSAM::INS;
        } else if (_v.type == Variant::DEL) {
            v.type = VariantFromSAM::DEL;
        }
        return v;
    }
    
    static inline std::vector<Variant> variantFromSAMVecToVariantVec(const std::vector<VariantFromSAM> &variants) {
        std::vector<Variant> newVariants;
        std::transform(variants.begin(), variants.end(), back_inserter(newVariants), variantFromSAMToVariant);
        return newVariants;
    }
    
    static inline std::vector<VariantFromSAM> variantVecToVariantFromSAMVec(const std::vector<Variant> &variants) {
        std::vector<VariantFromSAM> newVariants;
        std::transform(variants.begin(), variants.end(), back_inserter(newVariants), variantToVariantFromSAM);
        return newVariants;
    }
    
    static inline int getMinDistance(Variant &target, std::vector<Variant> &closeVariants) {
        int min = 10e5;
        for (std::vector<Variant>::iterator it = closeVariants.begin(); it != closeVariants.end(); it++) {
            int distance = std::abs(target.startInGenome - it->startInGenome);
            if (distance < min) {
                min = distance;
            }
        }
        return min;
    }
    
    static inline const Variant &findClosestVariant(const Variant &target, const std::vector<Variant> &closeVariants) {
        int min = 10e5;
        const Variant *closest;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin(); it != closeVariants.end(); it++) {
            int distance = std::abs(target.startInGenome - it->startInGenome);
            if (distance < min) {
                min = distance;
                closest = &(*it);
            }
        }
        return *closest;
    }
    
    static inline const Variant &findFurthestVariant(const Variant &target, const std::vector<Variant> &closeVariants) {
        int max = 0;
        const Variant *furthest;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin(); it != closeVariants.end(); it++) {
            int distance = std::abs(target.startInGenome - it->startInGenome);
            if (distance > max) {
                max = distance;
                furthest = &(*it);
            }
        }
        return *furthest;
    }
    
    static std::string getSymbols(std::vector<Variant> &variants) {
        if (variants.empty()) {
            return "-";
        } else if (variants.size() == 1) {
            return variants[0].getSymbol();
        }
        std::string out = variants[0].getSymbol();
        for (int i = 1;i < variants.size();i++) {
            out += "," + variants[i].getSymbol();
        }
        return out;
    }
    
    static inline std::pair<int, int> getRange(const Variant &target, const std::vector<Variant> &closeVariants) {
        int left = target.startInGenome, right = target.endInGenome;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin();it != closeVariants.end();it++) {
            if (left > it->startInGenome) {
                left = it->startInGenome;
            }
            if (right < it->endInGenome) {
                right = it->endInGenome;
            }
        }
        return std::make_pair(left, right);
    }
    
    static inline std::vector<Variant> excludeFurthestVariant(const Variant &target, const std::vector<Variant> &closeVariants) {
        if (closeVariants.size() <= 1) {
            return std::vector<Variant>();
        }
        const Variant &furthest = findFurthestVariant(target, closeVariants);
        std::vector<Variant> variants;
        for (std::vector<Variant>::const_iterator it = closeVariants.begin();it != closeVariants.end();it++) {
            if (furthest.startInGenome != it->startInGenome) {
                variants.push_back(*it);
            } else {
                LOG(logWARNING) << "exclude: " << *it << std::endl;
            }
        }
        return variants;
    }
    
    //
    // Utilities for candidate variants when checking pileup files
    //
    static inline bool isHeteroSNP(int ref, int obs, double confidenceProb){
        if(ref+obs > 0){
            return Utils::withinConfidenceInterval(ref+obs, obs, 0.50, confidenceProb);
        }else return false;
    }
    
    
    static inline void makeWindowsFromSortedCandidates(std::vector<CandidateVariant> &variants, std::vector<CandidateWindow> &windows,
    double log10_thres, int minDistToIndel, int maxWindowSize){
        windows.clear();
        std::vector<CandidateWindow>().swap(windows);
        windows.resize(variants.size());
        VariantUtils::fillHeteroSNPInfo(variants,log10_thres);
        VariantUtils::fillIndelCheckCover(variants, minDistToIndel);
        
        // init window
        for(int i = 0 ; i < variants.size(); i++){
            windows[i] = CandidateWindow(variants[i]);
        }
        // add close heteroSNP
        for(int i = 0 ; i < variants.size(); i++ ){
            if(variants[i].detail.isHeteroSNP) continue;
            int i_pos = variants[i].toPileUpStart(variants[i].type, variants[i].startInGenome);
            // down
            for(int j = i-1 ; 0 <= j ; j--){
                int j_pos = variants[j].toPileUpStart(variants[j].type, variants[j].startInGenome);
                if(variants[i].Chr != variants[j].Chr) break;
                if( abs(i_pos - j_pos) > maxWindowSize ) break;
                if(variants[j].detail.isHeteroSNP){
                    windows[i].addVariant(variants[j],true);
                }
            }            
            // up
            for(int j = i+1 ; j < variants.size(); j++){
                int j_pos = variants[j].toPileUpStart(variants[j].type, variants[j].startInGenome);
                if(variants[i].Chr != variants[j].Chr) break;
                if( abs(i_pos - j_pos) > maxWindowSize ) break;
                if(variants[j].detail.isHeteroSNP){
                    windows[i].addVariant(variants[j]);
                }
            }
        }
    }
    static inline void removeFalseHeteroSNP( std::vector<CandidateVariant> &candidates, double confidenceProb){
        int numSNP = VariantUtils::fillHeteroSNPInfo(candidates, confidenceProb);
        std::vector<CandidateVariant> ans(numSNP);
        int j = 0;
        for(int i = 0; i < candidates.size(); i++)if(candidates[i].detail.isHeteroSNP){
            ans[j] = candidates[i]; j++;
        }
        swap(ans,candidates);
    }
    
    static inline int fillHeteroSNPInfo( std::vector<CandidateVariant> &candidates, double confidenceProb){
        // check all variants and flag as heteroSNP or not.
        int ans = 0 ;
        for(int i = 0 ; i < candidates.size(); i++ ){
            int ref = candidates[i].detail.refNumN;
            int obs = candidates[i].detail.obsNumN;
            if( VariantUtils::isHeteroSNP( ref, obs, confidenceProb) ){
                candidates[i].detail.isHeteroSNP = true;
                ans ++;
            }else{
                candidates[i].detail.isHeteroSNP = false;
            }
        }
        return ans;
    }
    
    static inline void fillIndelCheckCover( std::vector<CandidateVariant> &candidates, int minDistToIndel ){
        for(int i = 0 ; i < candidates.size(); i++ ){
            int i_pos = candidates[i].startInGenome;
            candidates[i].detail.indelCoverCheck = false;
            // down
            for(int j = i-1 ; 0 <= j ; j--){
                if(candidates[i].Chr != candidates[j].Chr) break;
                int j_pos = candidates[j].startInGenome;
                if( abs(i_pos - j_pos) > minDistToIndel ) break;
                if(candidates[j].type == Variant::INS || candidates[j].type == Variant::DEL){
                    candidates[i].detail.indelCoverCheck = true;
                }
            }
            // up
            for(int j = i+1 ; j < candidates.size(); j++){
                if(candidates[i].Chr != candidates[j].Chr) break;
                int j_pos = candidates[j].startInGenome;
                if( abs(i_pos - j_pos) > minDistToIndel ) break;
                if(candidates[j].type == Variant::INS || candidates[j].type == Variant::DEL){
                    candidates[i].detail.indelCoverCheck = true;
                }
            }
        }
    }
    
    static inline void addToGlobalCandidates( std::vector<CandidateVariant> &all, std::vector<CandidateVariant> &local  ){
        if( local.size() > 1 ){
            for(int i = 0 ; i < local.size(); i++) local[i].detail.isTriallilic = true;
        }else {
            for(int i = 0 ; i < local.size(); i++) local[i].detail.isTriallilic = false;
        }
        for(int i = 0; i < local.size(); i++) all.push_back(local[i]);
    }
};

#endif /* defined(__haplotype_builder__VariantUtils__) */
