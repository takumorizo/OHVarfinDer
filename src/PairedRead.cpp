#include "PairedRead.h"

PairedRead::PairedRead(SamRead *read, SamRead *mateRead) {
    // first = read;
    // second = mateRead;
    SamRead*  plus = NULL;
    SamRead* minus = NULL;

    if(     read != NULL && !read->onReverseStrand )  plus=read;
    else if(read != NULL &&  read->onReverseStrand ) minus=read;

    if(     mateRead != NULL && !mateRead->onReverseStrand )  plus=mateRead;
    else if(mateRead != NULL &&  mateRead->onReverseStrand ) minus=mateRead;
    first = plus; second = minus;
}

std::string PairedRead::getSeqName() const {
    if(      first != NULL) return  first->seq_name;
    else if(second != NULL) return second->seq_name;
}

std::vector<IRead *> PairedRead::toIReadVector(std::vector<PairedRead> &reads) {
    std::vector<IRead *> ireads;
    ireads.reserve(reads.size());
    for (std::vector<PairedRead>::iterator it = reads.begin(); it != reads.end(); it++) {
        ireads.push_back(&(*it));
    }
    return ireads;
}

bool PairedRead::coverPosition(int pos) const {
    if (first  != NULL &&  first->cover(pos)) return true;
    if (second != NULL && second->cover(pos)) return true;
    return false;
}

bool PairedRead::hasVariant(const VariantFromSAM &variant) const {
    if (first  != NULL &&  first->hasVariantFromSam(variant)) return true;
    if (second != NULL && second->hasVariantFromSam(variant)) return true;
    return false;
}

double PairedRead::getMapQuality() const {
    bool isFirstOK  = (first  != NULL &&  !first->isUnmapped());
    bool isSecondOK = (second != NULL && !second->isUnmapped());
    if(       isFirstOK &&  isSecondOK)  return  first->mapQual * second->mapQual;
    else if(  isFirstOK && !isSecondOK)  return  first->mapQual;
    else if( !isFirstOK &&  isSecondOK)  return second->mapQual;
    else if( !isFirstOK && !isSecondOK)  return 0;
    else                                 return 0;
}

bool PairedRead::hasIndel() const {
    bool isFirstOK  = (first  != NULL &&  !first->isUnmapped());
    bool isSecondOK = (second != NULL && !second->isUnmapped());
    if(       isFirstOK &&  isSecondOK)  return  first->hasIndel || second->hasIndel;
    else if(  isFirstOK && !isSecondOK)  return  first->hasIndel;
    else if( !isFirstOK &&  isSecondOK)  return second->hasIndel;
    else if( !isFirstOK && !isSecondOK)  return false;
    else                                 return false;
}

bool PairedRead::hasSoftClip() const {
    bool isFirstOK  = (first  != NULL &&  !first->isUnmapped());
    bool isSecondOK = (second != NULL && !second->isUnmapped());
    if(       isFirstOK &&  isSecondOK)  return  first->hasSoftClip || second->hasSoftClip;
    else if(  isFirstOK && !isSecondOK)  return  first->hasSoftClip;
    else if( !isFirstOK &&  isSecondOK)  return second->hasSoftClip;
    else if( !isFirstOK && !isSecondOK)  return false;
    else                                 return false;
}

bool PairedRead::hasIndelPrimary() const {
    bool isFirstOK  = (first  != NULL &&  !first->isUnmapped());
    bool isSecondOK = (second != NULL && !second->isUnmapped());
    if(isFirstOK)  return   first->hasIndel;
    if(isSecondOK) return  second->hasIndel;
    return false;
}

bool PairedRead::hasSoftClipPrimary() const {
    bool isFirstOK  = (first  != NULL &&  !first->isUnmapped());
    bool isSecondOK = (second != NULL && !second->isUnmapped());
    if(isFirstOK)  return   first->hasSoftClip;
    if(isSecondOK) return  second->hasSoftClip;
    return false;
}

double PairedRead::getMapQualityPrimary() const {
    bool isFirstOK  = (first  != NULL &&  !first->isUnmapped());
    bool isSecondOK = (second != NULL && !second->isUnmapped());
    if(isFirstOK)  return   first->mapQual;
    if(isSecondOK) return  second->mapQual;
    return false;
}

SamRead *PairedRead::getPrimarySamRead(int pos) const {
    bool isFirstNonNull  = (first  != NULL);
    bool isSecondNonNull = (second != NULL);
    bool isFirstCover  = false;    if(isFirstNonNull)  isFirstCover  = first->cover(pos);
    bool isSecondCover = false;    if(isSecondNonNull) isSecondCover = second->cover(pos);

    if(        isFirstCover && isSecondCover){
        if(      first->hasFilter(64)){ return first;}
        else if(second->hasFilter(64)){ return second;}
    }else if(  isFirstCover && !isSecondCover){
        return first;
    }else if( !isFirstCover &&  isSecondCover){
        return second;
    }else if( !isFirstCover && !isSecondCover){
        return NULL;
    }else{
        return NULL;
    }
}

SamRead *PairedRead::getPrimarySamRead() const {
    bool isFirstOK  = (first  != NULL);
    bool isSecondOK = (second != NULL);
    if(isFirstOK)  return   first;
    if(isSecondOK) return  second;
    return NULL;
}