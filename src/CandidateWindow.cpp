/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <cstdlib>
#include <sstream>
#include <algorithm>
#include "CandidateWindow.h"
#include "Utils.h"

// ex: "1495:G=>C", "1495:G=>C,1410:+CGG", "1495:-AA"
std::vector<Variant> CandidateWindow::parseCloseVariants(std::string str) {
    std::vector<Variant> variants;
    if (str == "") {
        throw 1;
    }
    if (str == "-") {
        return variants;
    }
    std::vector<std::string> var;
    std::vector<std::string> variantsInStr = Utils::split(str, ",");
    variants.reserve(variantsInStr.size());
    for (int i = 0; i < variantsInStr.size(); i++) {
        var = Utils::split(variantsInStr[i], ":");
        Variant v(var[1], atoi(var[0].c_str()));
        variants.push_back(v);
    }
    return variants;
}

inline std::string CandidateWindow::getVariantSymbol(std::string ref, std::string obs) {
    if (ref == "-") {
        return ("+" + obs);
    } else if (obs == "-") {
        return ("-" + ref);
    } else {
        return (ref + "=>" + obs);
    }
}

CandidateWindow::CandidateWindow(std::string line) {
    try {
        std::stringstream linestream(line);
        CandidateWindow::Information &vi = info;
        std::string closeGermlineVariants;
        linestream >> vi.chr >> vi.start >> vi.end >> vi.ref >> vi.obs
        >> vi.ref_count_tumor >> vi.obs_count_tumor >> vi.ref_count_normal
        >> vi.obs_count_normal >> vi.missrate_tumor >> vi.strandrate_tumor
        >> vi.missrate_normal >> vi.strandrate_normal >> vi.ref_bq_tumor
        >> vi.obs_bq_tumor >> vi.ref_bq_normal >> vi.obs_bq_normal;
        if ((int)std::count(line.begin(), line.end(), '\t') == 20) {
            linestream >> vi.triallelic_site_check >> vi.indel_cover_check
                       >> vi.fisher_score >> closeGermlineVariants;
        } else {
            linestream >> vi.fisher_score >> closeGermlineVariants;
        }
        target = Variant(getVariantSymbol(vi.ref, vi.obs), vi.start);
        closeVariants = parseCloseVariants(closeGermlineVariants);
    } catch (...) {
        throw std::string("parsing candidate window failed: " + line);
    }
}


std::string CandidateWindow::toString(){
    std::string ans = "";
    ans += info.chr;                    ans += "\t";
    ans += Utils::toString(info.start); ans += "\t";
    ans += Utils::toString(info.end);   ans += "\t";
    ans += info.ref;                    ans += "\t";
    ans += info.obs;                    ans += "\t";
    ans += info.ref_count_tumor;        ans += "\t";
    ans += info.obs_count_tumor;        ans += "\t";
    ans += info.ref_count_normal;       ans += "\t";
    ans += info.obs_count_normal;       ans += "\t";
    ans += info.missrate_tumor;         ans += "\t";
    ans += info.strandrate_tumor;       ans += "\t";
    ans += info.missrate_normal;        ans += "\t";
    ans += info.strandrate_normal;      ans += "\t";
    ans += info.ref_bq_tumor;           ans += "\t";
    ans += info.obs_bq_tumor;           ans += "\t";
    ans += info.ref_bq_normal;          ans += "\t";
    ans += info.obs_bq_normal;          ans += "\t";
    ans += info.triallelic_site_check;  ans += "\t";
    ans += info.indel_cover_check;      ans += "\t";
    ans += info.fisher_score;           ans += "\t";
    
    for(int i = 0; i < closeVariants.size(); i++ ){
        ans += closeVariants[i].getSymbol();
        if(i != closeVariants.size()-1) ans += ",";
    }
    return ans;
}

CandidateWindow::CandidateWindow(const CandidateWindow &cand) {
    info = cand.info;
    target = cand.target;
    closeVariants = cand.closeVariants;
}

CandidateWindow &CandidateWindow::operator=(const CandidateWindow &cand) {
    info = cand.info;
    target = cand.target;
    closeVariants = cand.closeVariants;
    return *this;
}
/*
struct Information {
    std::string chr;
    int start, end;
    std::string ref, obs;
    std::string ref_count_tumor, obs_count_tumor;
    std::string missrate_tumor, strandrate_tumor;
    std::string ref_count_normal, obs_count_normal;
    std::string missrate_normal, strandrate_normal;
    std::string ref_bq_tumor, obs_bq_tumor, ref_bq_normal, obs_bq_normal;
    std::string triallelic_site_check, indel_cover_check;
    std::string fisher_score;
    //  std::string estimated_error_rate;
};
Information info;
Variant target;
std::vector<Variant> closeVariants;
*/

CandidateWindow::CandidateWindow(const CandidateVariant& candidate){
    info.chr   = candidate.Chr;
    info.start = candidate.startInGenome;
    info.end   = candidate.endInGenome;
    info.ref   = candidate.ref;
    info.obs   = candidate.obs;
    
    info.ref_count_tumor  = Utils::toString(candidate.detail.refNumT);
    info.obs_count_tumor  = Utils::toString(candidate.detail.obsNumT);
    info.ref_count_normal = Utils::toString(candidate.detail.refNumN);
    info.obs_count_normal = Utils::toString(candidate.detail.obsNumN);

    info.missrate_tumor    = Utils::toString(candidate.detail.rateT);
    info.missrate_normal   = Utils::toString(candidate.detail.rateN);
    info.strandrate_normal = Utils::toString(-1.0);
    info.strandrate_tumor  = Utils::toString(-1.0);
    info.ref_bq_tumor      = Utils::toString(-1.0);
    info.obs_bq_tumor      = Utils::toString(-1.0);
    info.ref_bq_normal     = Utils::toString(-1.0);
    info.obs_bq_normal     = Utils::toString(-1.0);
    
    if(candidate.detail.isTriallilic) info.triallelic_site_check = "triallelic";
    else                              info.triallelic_site_check = "ok";
    
    if(candidate.detail.indelCoverCheck) info.indel_cover_check = "indel-cover";
    else                                 info.indel_cover_check = "ok";
    
    info.fisher_score = "-";
    target = candidate;
}

void CandidateWindow::addVariant(CandidateVariant candidate,bool isFront){
    if(isFront) closeVariants.insert( closeVariants.begin(),((Variant)candidate));
    else closeVariants.push_back((Variant)candidate);
}

