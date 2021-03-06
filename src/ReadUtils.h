#ifndef __haplotype_builder__ReadUtils__
#define __haplotype_builder__ReadUtils__

#include <vector>
#include <string>
#include <algorithm>

#include "Variant.h"
#include "IRead.h"
class ReadUtils {
public:
    static inline IRead *findBySeqName(const std::vector<IRead *> &reads, std::string seqName) {
        for (std::vector<IRead *>::const_iterator it = reads.begin(); it != reads.end(); it++) {
            if ((*it)->getSeqName() == seqName)
                return *it;
        }
        throw std::string("not_found");
    }
    
    static inline std::vector<IRead *> coverPosition(const std::vector<IRead *> &reads, int pos) {
        std::vector<IRead *> vec;
        for (std::vector<IRead *>::const_iterator it = reads.begin(); it != reads.end(); it++) {
            if ((*it)->coverPosition(pos)) {
                vec.push_back(*it);
            }
        }
        return vec;
    }
    
    static inline std::vector<IRead *> hasVariant(const std::vector<IRead *> &reads, VariantFromSAM &v) {
        std::vector<IRead *> vec;
        for (std::vector<IRead *>::const_iterator it = reads.begin(); it != reads.end(); it++) {
            if ((*it)->hasVariant(v)) {
                vec.push_back(*it);
            }
        }
        return vec;
    }
    
    static inline std::vector<IRead *> filterByMapQuality(const std::vector<IRead *> &reads, double threshold) {
        std::vector<IRead *> vec;
        for (std::vector<IRead *>::const_iterator it = reads.begin(); it != reads.end(); it++) {
            if ((*it)->getMapQuality() >= threshold) {
                vec.push_back(*it);
            }
        }
        return vec;
    }
    
    static inline bool coverVariants(const IRead* read, const std::vector<Variant>& variants ){
        bool isCover = false;
        for(int i = 0; i < variants.size(); i++){
            int pos = variants[i].startInGenome;
            isCover = isCover || read->coverPosition(pos);
            if(isCover) break;
        }
        return isCover;
    }
    
};

#endif /* defined(__haplotype_builder__ReadUtils__) */
