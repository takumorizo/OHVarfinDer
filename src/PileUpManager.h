//
//  PileUpManager.hpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/05/24.
//  Copyright © 2016年 森山卓也. All rights reserved.
//

#ifndef PileUpManager_h
#define PileUpManager_h

#include <vector>
#include "CandidateWindow.h"
#include "Parameters.h"
#include "PileUpUtils.h"
#include "PileUp.h"
#include "CandidateVariant.h"
#include "PileUpChecker.h"

class PileUpManager{
public:
    PileUpManager();
    std::vector<CandidateWindow> searchCandidateWindow(Parameters param);
private:
	bool isHeteroSNP(const std::vector<std::string> &cols, PileUp &pileN, 
					 std::vector<CandidateVariant> &tempCandidates, 
					 const Parameters &param, PileUpChecker &pileChecker);
	bool isTriallelic(const  std::vector<std::string> &cols, PileUp &pileT,
					  std::vector<CandidateVariant> &tempCandidates, 
					  const Parameters &param, PileUpChecker &pileChecker);
	bool isSomaticCandidate(const  std::vector<std::string> &cols, PileUp &pileN, PileUp &pileT,
                            std::vector<CandidateVariant> &tempCandidates, 
                            const Parameters &param, PileUpChecker &pileCheckerN, PileUpChecker &pileCheckerT);

};


#endif /* PileUpManager_h */
