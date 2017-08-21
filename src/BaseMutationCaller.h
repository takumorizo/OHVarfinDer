//
//  BaseMutationCaller.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/04/06.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef __OHVarFinder__BaseMutationCaller__
#define __OHVarFinder__BaseMutationCaller__

#include "faidx.h"
#include "BamReader.h"
#include "Parameters.h"
#include "MutationCallResult.h"
#include "CandidateWindow.h"
#include "HaplotypeBuilder.h"
#include <sstream>

class BaseMutationCaller {
public:
//    HaplotypeBuilder hapBuilder;
//    BamReader &tumorBamReader;
//    BamReader &normalBamReader;
//    Parameters &params;
//    const faidx_t *fai;
    virtual ~BaseMutationCaller(){};
    BaseMutationCaller(){};
//    BaseMutationCaller(BamReader &_tumorBamReader, BamReader &_normalBamReader, Parameters &_params,const faidx_t *_fai) : params(_params), tumorBamReader(_tumorBamReader), normalBamReader(_normalBamReader), hapBuilder(_fai), fai(_fai) {
//        hapBuilder.windowSize = 2 * (params.maxInsertSize + params.maxReadLength);
//    }
    virtual MutationCallResult call(const CandidateWindow &candidateWindow) {
        MutationCallResult result(candidateWindow) ; return result;
    }
};

#endif /* defined(__OHVarFinder__BaseMutationCaller__) */
