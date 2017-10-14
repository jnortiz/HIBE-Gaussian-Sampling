#include "KnuthYaoSamplerPARI.h"

int main(){
    pari_init(10000000, 2);
    globalvars* g = new globalvars;
    parameters* params = new parameters;
    params->s = 4;
    params->sigma = 6;
    g->pPack = genProbabilityMatrix(params, "0.00");
    printf("Probobility matrix generated\n");
    int sampledvalue = SampleKnuthYao(0, params, g);
    printf("%d\n", sampledvalue);
    pari_close();
}
