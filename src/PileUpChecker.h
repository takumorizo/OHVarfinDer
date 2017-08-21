//
//  PileUpChecker.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/14.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef __OHVarFinder__PileUpChecker__
#define __OHVarFinder__PileUpChecker__

#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#include <ctype.h>
#include "Utils.h"
#include <iostream>
#include "Parameters.h"
#include "PileUpUnit.h"
#include "PileUp.h"
#include "CandidateVariant.h"

class PileUpChecker {
public:
    Parameters::PileUpParameters param;
    PileUpChecker(Parameters::PileUpParameters param);

    void setCandidateMutations(std::vector<CandidateVariant> &candidates, PileUp &pileUp,
                               const std::string &Chr,    const std::string &position,
                               const std::string &ref,    const std::string &depth,
                               const std::string &bases,  const std::string &qualities);
    void filterCandidates(std::vector<CandidateVariant> &candidates,const PileUp &pileUp, bool isTumor);
    
};

#endif /* defined(__OHVarFinder__PileUpChecker__) */


