#ifndef ReadTest_Utils_h
#define ReadTest_Utils_h

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <ctype.h>
#include <sstream>

class Utils {
public:
    static inline std::vector<std::string> split(std::string s, std::string delimiter) {
        std::vector<std::string> vec;
        size_t pos = 0;
        std::string token;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            vec.push_back(token);
            s.erase(0, pos + delimiter.length());
        }
        if (s.size() != 0) {
            vec.push_back(s);
        }
        return vec;
    }
    
    static inline bool nearlyEqual(double a, double b) {
        return std::abs(a - b) < 0.01;
    }
    
    static inline int countChar(const std::string &s, char c){
        int ans = 0;
        for(int i = 0; i < s.length(); i++){
            if(s[i] == c) ans++;
        }
        return ans;
    };
    
    static inline void toUpperString( std::string &s ){
        for(int i = 0 ; i < s.length(); i++) s[i] = toupper(s[i]);
    }
    static inline void toLowerString( std::string &s ){
        for(int i = 0 ; i < s.length(); i++ ) s[i] = tolower(s[i]);
    }
    
    template<class T> static inline std::string toString(T x) {std::ostringstream sout;sout<<x;return sout.str();}
    
#include "float.h"
#include "math.h"
    // naively counts the number of supporting reads
    static inline std::vector<double> votes(const std::vector<std::vector<double> > &liks) {
        if (liks.size() == 0) {
            throw std::string("liks.size() == 0");
        }
        std::vector<double> vs = std::vector<double>(liks[0].size());
        for (int i = 0; i < liks.size(); i++) {
            double norm = 0.0;
            for (int j = 0; j < liks[i].size(); j++) {
                norm += exp( liks[i][j] );
            }
            for (int j = 0; j < liks[i].size(); j++) {
                vs[j] += exp( liks[i][j] - log(norm) );
            }
        }
        return vs;
    }
    
    // naively counts the number of supporting reads, in the index vector.
    static inline std::vector<double> votes(const std::vector<std::vector<double> > &liks, std::vector<int> idx){
        if (liks.size() == 0) {
            throw std::string("liks.size() == 0");
        }
        std::vector<double > vs = std::vector<double>(idx.size());
        for(int i = 0; i < liks.size(); i++){
            double norm = 0.0;
            for(int j = 0; j < idx.size(); j++){
                norm += exp(liks[i][idx[j]]);
            }
            for(int j = 0; j < idx.size(); j++){
                vs[j] += exp( liks[i][idx[j]] - log(norm) );
            }
        }
        return vs;
    }
    
    // naively counts the number of supporting reads, in the index vector.
    static inline std::vector<double> votes(const std::vector<std::vector<double> > &liks,std::vector<std::vector<int> > types, std::vector<int> idx){
        if (liks.size() == 0) {
            throw std::string("liks.size() == 0");
        }
        std::vector<double > vs = std::vector<double>(idx.size());
        for(int i = 0; i < liks.size(); i++)if(types[i][idx[0]] != 0 ){
            double norm = 0.0;
            for(int j = 0; j < idx.size(); j++){
                norm += exp(liks[i][idx[j]]);
            }
            for(int j = 0; j < idx.size(); j++){
                vs[j] += exp( liks[i][idx[j]] - log(norm) );
            }
        }
        return vs;
    }
    
    static inline double log10_Prob_Binomial(int n, int k, double p){
        if(n == 0) throw std::string("unexpected : n == 0");
        if(p< 0.0 || p > 1.0) throw std::string("unexpected : p < 0 or p > 1");
        //log10_n_C_k
        double log10_n_C_k = 0.0;
        for(int i = 1 ; i <= n     ; i++) log10_n_C_k += std::log10(1.0*i);
        for(int i = 1 ; i <= (n-k) ; i++) log10_n_C_k -= std::log10(1.0*i);
        for(int i = 1 ; i <= k     ; i++) log10_n_C_k -= std::log10(1.0*i);
        log10_n_C_k +=     k*std::log10(p);
        log10_n_C_k += (n-k)*std::log10((1.0-p));
        return log10_n_C_k;
    }
    
    static inline bool withinConfidenceInterval(int n, int k, double p, double confidenceProb){
        int mid = (int)(n * p);
        double interval_prob = pow(10, log10_Prob_Binomial(n, mid, p));
        std::map<int, bool > confidentInt;
        confidentInt[mid] = true;
        int i = 1;
        while(interval_prob < confidenceProb){
            interval_prob += pow(10, log10_Prob_Binomial(n, mid+i, p));
            confidentInt[mid+i] = true;
            if(mid -i >= 0){
                interval_prob += pow(10,log10_Prob_Binomial(n, mid-i, p));
                confidentInt[mid-i] = true;
            }
            i++;
        }
        return confidentInt.count(k);
    }

};

#endif
