/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#ifndef __BamReader__
#define __BamReader__

#include <vector>
#include <string>

#include "MyBam.hpp"
#include "SamRead.h"
#include "PairedRead.h"
#include "Parameters.h"

class BamReader {
    // readBuffer is responsible for the memory management of Read *
    std::vector<SamRead *> readBuffer;
    const MyBam *bam;
    const faidx_t *fai;
    const int fFlag;
    const int FFlag;
public:
    int maxInsertSize, maxReads;
    BamReader(const MyBam *bam, const faidx_t *, Parameters param,
              int _maxReads = 20000);
    std::vector<IRead *> getPairedReads(std::string chr, int targetPos,
                                        int left, int right, int maxInsertSize = 1000);
    std::vector<IRead *> getSingleReads(std::string chr, int targetPos,
                                        int left, int right);
    ~BamReader();
};

#endif /* defined(__BamReader__) */
