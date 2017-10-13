#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <pari/pari.h>
#include <time.h>
#define PARI_OLD_NAMES
#include <vector>
#include <sys/time.h>

#define precision 6 // This means 256-bit precision

struct parameters{
    int s;
    int sigma;
};

struct ProbMatrixPack{
    GEN P;
    std::vector<int> startPos;
    bool isInitialized = false;
} *pPackglobal;

struct globalvars{
    ProbMatrixPack* pPack;
};

GEN getGuassProbability(GEN point, GEN center, parameters* params){
    int sigma = params->sigma;
    GEN twopi = mulir(stoi(2), mppi(precision));
    GEN s = mulir(stoi(sigma), sqrtr(twopi));
    GEN sinv = invr(s);
    if(gcmp(point, strtor("0.00", precision)) == 0)
        return sinv;
    else{
        return gmul(sinv, gexp( gdiv(gneg(gpow(gdiv(gsub(point, center), stoi(sigma)), stoi(2), precision)), strtor("2.00", precision)), precision));
    }
}

ProbMatrixPack* genProbabilityMatrix(parameters* params, char* c){
    int sigma = params->sigma;
    int tailprune = params->s;
    GEN center;
    center = strtor(c, precision);
    GEN tempP, beginAddressP, ProbofPoints;
    int bounds = tailprune*sigma;
    ProbofPoints = cgetg(bounds+2, t_REAL);
    int bitprecision = 64*(precision-2);
    tempP = cgetg(bitprecision+1, t_VEC);
    for(int i=1; i<=bitprecision; i++){
        GEN temp = cgetg(bounds+2, t_INT);
        gel(tempP, i) = temp;
    }
    
    for(int x = bounds; x > 0; x--){
        gel(ProbofPoints, bounds+1-x) = getGuassProbability(gadd(center, stoi(x)), center, params);
    }
    gel(ProbofPoints, bounds+1) = gdiv(getGuassProbability(gadd(center, stoi(0)), center, params), stoi(2));
    
    int i = -1;
    for(int j=0; j<bitprecision; j++){
        GEN temppow = gpow(strtor("2.00", precision), stoi(i), precision);
        i--;
        for(int x = bounds; x >= 0; x--){
            gel(gel(tempP, j+1), bounds+1-x) = stoi(0);
            if(gcmp(gel(ProbofPoints, bounds+1-x), temppow) >= 0){
                gel(gel(tempP, j+1), bounds+1-x) = stoi(1);
                gel(ProbofPoints, bounds+1-x) = gsub(gel(ProbofPoints, bounds+1-x), temppow);
            }
        }
    }
    std::vector<int> beginPos;
    
    for(int x = bounds; x >= 0; x--){
        for(int j=0; j<bitprecision; j++){
            if(j == bitprecision-2){
                beginPos.push_back(j);
                break;
            }
            
            if(gcmp(gel(gel(tempP, j+1), bounds+1-x), stoi(1))==0){
                beginPos.push_back(j);
                break;
            }
        }
    }
    ProbMatrixPack* pPack = new ProbMatrixPack;
    GEN P;
    P = tempP;
    pPack->P = P;
    pPack->startPos = beginPos;
    pPack->isInitialized = true;
    pPackglobal = pPack;
    return pPack;
    
}

ProbMatrixPack* genProbabilityMatrix(parameters* params, int c){
    int sigma = params->sigma;
    int tailprune = params->s;
    GEN center;
    center = stoi(c);
    GEN tempP, beginAddressP, ProbofPoints;
    
    int bounds = tailprune*sigma;
    ProbofPoints = cgetg(bounds+2, t_REAL);
    int bitprecision = 64*(precision-2);
    tempP = cgetg(bitprecision+1, t_VEC);
    for(int i=1; i<=bitprecision; i++){
        GEN temp = cgetg(bounds+2, t_INT);
        gel(tempP, i) = temp;
    }
    for(int x = bounds; x > 0; x--){
        gel(ProbofPoints, bounds+1-x) = getGuassProbability(gadd(center, stoi(x)), center, params);
    }
    gel(ProbofPoints, bounds+1) = gdiv(getGuassProbability(gadd(center, stoi(0)), center, params), stoi(2));

    
    int i = -1;
    for(int j=0; j<bitprecision; j++){
        GEN temppow = gpow(strtor("2.00", precision), stoi(i), precision);
        i--;
        for(int x = bounds; x >= 0; x--){
            gel(gel(tempP, j+1), bounds+1-x) = stoi(0);
            if(gcmp(gel(ProbofPoints, bounds+1-x), temppow) >= 0){
                gel(gel(tempP, j+1), bounds+1-x) = stoi(1);
                gel(ProbofPoints, bounds+1-x) = gsub(gel(ProbofPoints, bounds+1-x), temppow);
            }
        }
    }
    
    std::vector<int> beginPos;
    
    for(int x = bounds; x >= 0; x--){
        for(int j=0; j<bitprecision; j++){
            if(j == bitprecision-2){
                beginPos.push_back(j);
                break;
            }
            
            if(gcmp(gel(gel(tempP, j+1), bounds+1-x), stoi(1))==0){
                beginPos.push_back(j);
                break;
            }
        }
    }
    ProbMatrixPack* pPack = new ProbMatrixPack;
    GEN P;
    P = tempP;
    pPack->P = P;
    pPack->startPos = beginPos;
    pPack->isInitialized = true;
    pPackglobal = pPack;
    return pPack;
    
}

int SampleKnuthYao(int c, parameters* params, globalvars* g){
    int sigma = params->sigma;
    int tailprune = params->s;
    GEN center;
    center = stoi(c);
    
    int bounds, col, d, invsample, pRows, pCols, s, flag, enable, hit;
    unsigned long r;
    bounds = tailprune*sigma;
    d = 0;
    hit = 0;
    invsample = bounds+1;
    
    GEN P = g->pPack->P;
    std::vector<int> beginPos = g->pPack->startPos;
    int bitprecision = 64*(precision-2);
    pRows = lg(P)-1;
    pCols = bitprecision;
    
    flag = 1-2*(rand()%2);
    
    int randomBits[pRows];
    int length = sizeof(unsigned long)*8;
    
    
    for(int i=0; i<pRows; i++){
        randomBits[i] = rand()%2;
    }
    
    s = 0;
    enable = 0;
    for(int row = 0; row<pRows; row++){
        if(enable==1)
            break;
        d = 2*d + randomBits[row];
        for(int col = beginPos[row]; col < pCols; col++) {
            d = d - itos(gel(gel(P, col+1), row+1));
            if(d==-1){
                hit = 1;
                s = col;
                enable = 1;
                break;
            }
        }
        
    }
    return (s*flag)+itos(center);
}
