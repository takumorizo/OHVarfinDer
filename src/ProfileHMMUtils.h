#ifndef __haplotype_builder__ProfileHMMUtils__
#define __haplotype_builder__ProfileHMMUtils__

#include "float.h"

#include <typeinfo>
#include <vector>
#include <map>

#include "Alignment.h"
#include "PairedRead.h"
#include "SingleRead.h"
#include "Haplotype.hpp"
#include "ProfileHMM.h"
#include "Parameters.h"
#include <utility>
#include "log.h"
#include <stdexcept>

class ProfileHMMUtils {
public:
    static inline double calc(const std::vector<Variant> &allVariants,
                              const Haplotype &hap, const PairedRead &read) {
        double output = 0.0;
        bool covered1 = false, covered2 = false;
        for (int i = 0; i < allVariants.size(); i++) {
            if (read.first != NULL  && read.first->cover(allVariants[i].startInGenome)) {
                covered1 = true;
            }
            if (read.second != NULL && read.second->cover(allVariants[i].startInGenome)) {
                covered2 = true;
            }
            if (covered1 && covered2) {
                break;
            }
        }
        if (covered1) {
            try {
                output += calc(hap, read.first);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "plus strand read mapping fail; continue;" << std::endl;
            }
        }
        if (covered2) {
            try {
                output += calc(hap, read.second);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "minus strand read mapping fail; continue;" << std::endl;
            }
        }
        return output;
    }

    static inline double calc(const std::vector<Variant> &allVariants,
                              const Haplotype &hap, const SingleRead &read) {
        bool covered = false;
        for (int i = 0; i < allVariants.size(); i++) {
            if (read.read->cover(allVariants[i].startInGenome)) {
                covered = true;
                break;
            }
        }
        if (covered) {
            return calc(hap, read.read);
        } else {
            return 0.0;
        }
    }
    
    static inline double calc(const std::vector<Variant> &allVariants,
                              const Haplotype &hap, const IRead *read) {
        if (typeid(*read) == typeid(PairedRead)) {
            return calc(allVariants, hap, *dynamic_cast<const PairedRead *>(read));
        } else if (typeid(*read) == typeid(SingleRead)) {
            return calc(allVariants, hap, *dynamic_cast<const SingleRead *>(read));
        } else {
            throw std::string("ProfileHMMUtils::calc not_implemented");
        }
    }
    
    static inline double calc(const Haplotype &hap, const SamRead *sam) {
        try{
            ProfileHMM pHMM0 = genProfileHMM(hap, sam);
            Alignment a0 = pHMM0.viterbi();
            return a0.likelihood;
        } catch (...) {
            throw std::string("bug in profileHMM");
        }
        /* TODO
        int th1 = -100, th2 = -200;
        if (a0.likelihood > th1) {
           // a0.print();
            return a0.likelihood;
        } else {
            if (sam->hasSoftClip) {
                std::string newSeq;
                std::vector<double> newBQ;
                int newPos;
                sam->getReadWithoutSoftClipping(&newSeq, &newBQ, &newPos);
                ProfileHMM pHMM1 = ProfileHMM(hap, newSeq, newBQ, newPos,
                                              sam->mapQual,
                                              true,
                                              sam->onReverseStrand, sam->seq_name);
                Alignment a1 = pHMM1.viterbi();
                if (a1.likelihood > th2) {
                    if (a1.likelihood < th1) {
                        a1.print();
                    }
                    LOG(logDEBUG) << "aligned without soft-clip (" << sam->seq_name << "): " << a1.likelihood << std::endl;
                    return a1.likelihood;
                } else {
                    a0.print();
                    a1.print();
                    throw std::string("did_not_map_profileHMM2");
                }
            } else if (a0.likelihood > th2) {
                return a0.likelihood;
            } else {
                a0.print();
                throw std::string("did_not_map_profileHMM1");
            }
        }
         */
    }
    
    static std::vector<std::vector<double> > calcAllLikelihoods(const std::vector<Haplotype> &haps,
                                                                const std::vector<IRead *> &reads) {
        std::vector<std::vector<double> > liks;
        liks.reserve(reads.size());
        std::vector<Variant> allVariants = Haplotype::getAllVariants(haps);
    
        for (int i = 0; i < reads.size(); i++) {
            try {
                std::vector<double> tmpLiks(haps.size());
                for (int j = 0; j < haps.size(); j++) {
                    tmpLiks[j] = calc(allVariants, haps[j], reads[i]);
                }
                liks.push_back(tmpLiks);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }

    static inline int dot( const std::vector<int> &array, const std::vector<int> &indexes){
        int ans = 0;
        for(int i = 0; i < std::max( array.size(), indexes.size() ); i++ ){
            ans += array[i] * indexes[i];
        }
        return ans;
    }
    
    static inline std::vector<std::vector<double> > calcAllLikelihoodsDiff(const std::vector<std::pair<Haplotype, Haplotype> > &hapPairs,
                                                                           const std::vector<IRead *> &reads,
                                                                           const std::vector<std::vector<int> > &types,
                                                                           const std::map<std::pair<int,int>, std::pair<int,int> > &hapIndexes,
                                                                           const std::vector<std::vector<double> > &originalLiks,
                                                                           const std::vector<std::vector<int> > &typesOld){
        std::vector<std::vector<double> > liks;
        std::map<std::pair<int,int>, double > tmpAnsLikDic;

        for (int i = 0; i < reads.size(); i++){
            if( dot(types[i], typesOld[i]) != 0 ){ 
                liks.push_back(originalLiks[i]);
            }else if( dot(types[i], typesOld[i]) == 0 ){
                try {
                    if (typeid(*(reads[i])) == typeid(PairedRead)) {
                        liks.push_back(calcEfficient(hapPairs, *dynamic_cast<const PairedRead *>(reads[i]), types[i], hapIndexes, tmpAnsLikDic));
                    } else if (typeid(*(reads[i])) == typeid(SingleRead)) {
                        liks.push_back(calcEfficient(hapPairs, *dynamic_cast<const SingleRead *>(reads[i]), types[i], hapIndexes, tmpAnsLikDic));
                    }
                } catch (std::string &s) {
                    LOG(logDEBUG) << s << std::endl;
                    LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
                }
            }
        }
        return liks;
    }    

    // debug for single end read
    static inline std::vector<std::vector<double> > calcAllLikelihoods(const std::vector<std::pair<Haplotype, Haplotype> > &hapPairs,
                                                                       const std::vector<IRead *> &reads,
                                                                       const std::vector<std::vector<int> > &types,
                                                                       const std::map<std::pair<int,int>, std::pair<int,int> > &hapIndexes){
        std::vector<std::vector<double> > liks;
        std::map<std::pair<int,int>, double > tmpAnsLikDic;
        for (int i = 0; i < reads.size(); i++){
            try {
                if (typeid(*(reads[i])) == typeid(PairedRead)) {
                    liks.push_back(calcEfficient(hapPairs, *dynamic_cast<const PairedRead *>(reads[i]), types[i], hapIndexes, tmpAnsLikDic));
                } else if (typeid(*(reads[i])) == typeid(SingleRead)) {
                    liks.push_back(calcEfficient(hapPairs, *dynamic_cast<const SingleRead *>(reads[i]), types[i], hapIndexes, tmpAnsLikDic));
                }
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }

    static inline std::vector<double> calcEfficient(const std::vector<std::pair<Haplotype, Haplotype> > &hapPairs,
                                                    const PairedRead& read,
                                                    const std::vector<int> &types,
                                                    const std::map<std::pair<int,int>, std::pair<int,int> > &hapIndexes,
                                                    std::map<std::pair<int,int>, double > &ansLikDic){
        SamRead *rp = read.first;
        SamRead *rm = read.second;
        // LOG(logINFO) << "calcEfficient pair start " << std::endl;
        // std::map<std::pair<int,int>, double > ansLikDic;
        ansLikDic.clear();
        std::vector<double> tmpLiks(hapPairs.size(), 0);
        for (int j = 0; j < hapPairs.size(); j++) if(types[j] != 0) for(int p = 0; p < 2; p++) {
            std::pair<int,int> idx = std::make_pair(j,p);
            std::pair<int,int> actualIdx = idx;
            {
                std::map<std::pair<int,int>, std::pair<int,int> >::const_iterator it;
                it = hapIndexes.find(idx);
                if(it != hapIndexes.end()) {
                    actualIdx = it->second;
                }
            }
            double tmp = 0.0;
            if(       ansLikDic.count(actualIdx) == 0 ){// unique or first haplotype pattern
                if(      p == 0 && rp != NULL) tmp = calc(hapPairs[j].first,  rp); 
                else if( p == 1 && rm != NULL) tmp = calc(hapPairs[j].second, rm);
            }else if( ansLikDic.count(actualIdx) >  0 ){// non unique or previously existing haplotype pattern
                tmp = ansLikDic[actualIdx];
            }
            tmpLiks[j] += tmp; ansLikDic[actualIdx] = tmp; ansLikDic[idx] = tmp;
        }
        return tmpLiks;
    }
    // TODO : debug for single end read is not done.
    static inline std::vector<double> calcEfficient(const std::vector<std::pair<Haplotype, Haplotype> > &hapPairs,
                                                    const SingleRead& read,
                                                    const std::vector<int> &types,
                                                    const std::map<std::pair<int,int>, std::pair<int,int> > &hapIndexes,
                                                    std::map<std::pair<int,int>, double > &ansLikDic){
        SamRead *r = read.read;
        int strand = 0;
        if(r->onReverseStrand == false) strand = 1;
        // LOG(logINFO) << "calcEfficient single start " << std::endl;
        // std::map<std::pair<int,int>, double > ansLikDic;
        ansLikDic.clear();
        std::vector<double> tmpLiks(hapPairs.size(), 0);
        for (int j = 0; j < hapPairs.size(); j++) if(types[j] != 0){
            std::pair<int,int> idx = std::make_pair(j, strand);
            std::pair<int,int> actualIdx = idx;
            {
                std::map<std::pair<int,int>, std::pair<int,int> >::const_iterator it;
                it = hapIndexes.find(idx);
                if(it != hapIndexes.end()) {
                    actualIdx = it->second;
                }            
            }
            double tmp = 0.0;
            if(       ansLikDic.count(actualIdx) == 0 ){// unique or first haplotype pattern
                if(      strand == 0) tmp = calc(hapPairs[j].first,  r); 
                else if( strand == 1) tmp = calc(hapPairs[j].second, r);
            }else if( ansLikDic.count(actualIdx) >  0  ){// non unique or previously existing haplotype pattern
                tmp = ansLikDic[actualIdx];
            }
            tmpLiks[j] += tmp; ansLikDic[actualIdx] = tmp; ansLikDic[idx] = tmp;
        }
        return tmpLiks;
    }


    static inline std::vector<std::vector<double> > calcAllLikelihoods(const std::vector<std::pair<Haplotype, Haplotype> > &hapPairs,
                                                                       const std::vector<IRead *> &reads,
                                                                       const std::vector<std::vector<int> > &types){
//      std::pair<Haplotype, Haplotype> ; first is plus strand, second is minus strand
        std::vector<std::vector<double> > liks;
        for (int i = 0; i < reads.size(); i++) {
            try {
                std::vector<double> tmpLiks(hapPairs.size());
                for (int j = 0; j < hapPairs.size(); j++) {
                    if(     types[i][j] != 0) tmpLiks[j] = types[i][j]*calc(hapPairs[j], reads[i]);
                    else if(types[i][j] == 0) tmpLiks[j] = 0.0;
                }
                liks.push_back(tmpLiks);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }

    static inline double calc(const std::pair<Haplotype, Haplotype> &hapPair,
                              const PairedRead& read){
        double ans = 0.0;
        SamRead *rp = read.first;
        SamRead *rm = read.second;
        if( rp != NULL && rp->onReverseStrand ) std::swap(rp, rm);
        if( rp != NULL) ans += calc(hapPair.first, rp);
        if( rm != NULL) ans += calc(hapPair.second, rm);
        return ans;
    }
    //TODO : debug not done for single end read
    static inline double calc(const std::pair<Haplotype, Haplotype> &hapPair,
                              const SingleRead& read){
        double ans = 0.0;
        SamRead *r = read.read;
        if( r != NULL){
            if( r->onReverseStrand == false ) ans += calc(hapPair.first, r);
            else                              ans += calc(hapPair.second, r);
        }
        return ans;
    }
    //TODO : debug not done for single end read
    static inline double calc(const std::pair<Haplotype, Haplotype> &hapPair,
                              const IRead* read){
        if (typeid(*read) == typeid(PairedRead)) {
            return calc(hapPair, *dynamic_cast<const PairedRead *>(read));
        } else if (typeid(*read) == typeid(SingleRead)) {
            return calc(hapPair, *dynamic_cast<const SingleRead *>(read));
            
        } else {
            throw std::string("ProfileHMMUtils::calc not_implemented");
        }
    }
    
    
    
    
    static inline ProfileHMM genProfileHMM(const Haplotype &hap, const SamRead *sam) {
        return ProfileHMM(hap, sam->seq, sam->qual, sam->leftMostPos, sam->mapQual, sam->hasIndel || sam->hasSoftClip || sam->hasHardClip, sam->onReverseStrand, sam->seq_name);
    }
    
    static std::vector<std::vector<double> > calcLikelihoods(const std::vector<Haplotype> &haps,
                                                             const std::vector<IRead *> &reads) {
        std::vector<std::vector<double> > liks = calcAllLikelihoods(haps, reads);
        return liks;
    }

    static std::vector<std::vector<double> > calcBasicLikelihoods(const std::vector<Haplotype> &haps,
                                                                  const std::vector<IRead *> &reads, int pos) {
        std::vector<std::vector<double> > liks;
        liks.reserve(reads.size());

        for (int i = 0; i < reads.size(); i++) {
            try {
                std::vector<double> tmpLiks(haps.size());
                for (int j = 0; j < haps.size(); j++) {
                    tmpLiks[j] = calc(haps[j], reads[i]->getPrimarySamRead(pos));
                }
                liks.push_back(tmpLiks);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }

    static std::vector<std::vector<double> > calcBasicLikelihoods(const std::vector<Haplotype> &haps,
                                                                  const std::vector<IRead *> &reads) {
        std::vector<std::vector<double> > liks;
        liks.reserve(reads.size());

        for (int i = 0; i < reads.size(); i++) {
            try {
                std::vector<double> tmpLiks(haps.size());
                for (int j = 0; j < haps.size(); j++) {
                    tmpLiks[j] = calc(haps[j], reads[i]->getPrimarySamRead());
                }
                liks.push_back(tmpLiks);
            } catch (std::string &s) {
                LOG(logDEBUG) << s << std::endl;
                LOG(logWARNING) << "calc ProfileHMM fail. skip this read " << reads[i]->getSeqName() << std::endl;
            }
        }
        return liks;
    }

};



#endif /* defined(__haplotype_builder__ProfileHMMUtils__) */
