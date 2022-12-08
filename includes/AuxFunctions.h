#include <cstdio>
#include "TTree.h"

#define PBARSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBARW 60

void PrintProgress(float percentage) {
    int val = (int)(percentage*100);
    int lpad =(int)(percentage*PBARW);
    int rpad = PBARW-lpad;
    printf("\r%3d%% [%.*s%*s]",val,lpad,PBARSTR,rpad,"");
    fflush(stdout);
}

uint64_t SortEntries(uint64_t& firstEvent, uint64_t& maxEvents, TTree* h101) {
	firstEvent = std::min(firstEvent, (uint64_t)h101->GetEntries());
    uint64_t n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
	return n;
}

template<class T>
void ReleaseMalloc(T x) {
	free(x);
}
template<class T, class... Args>
void ReleaseMalloc(T x, Args... args) {
	free(x); ReleaseMalloc(args...);
}
