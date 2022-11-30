#include <cstdio>

#define PBARSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBARW 60

void PrintProgress(float percentage) {
    int val = (int)(percentage*100);
    int lpad =(int)(percentage*PBARW);
    int rpad = PBARW-lpad;
    printf("\r%3d%% [%.*s%*s]",val,lpad,PBARSTR,rpad,"");
    fflush(stdout);
}

template<class T>
void ReleaseMalloc(T x) {
	free(x);
}
template<class T, class... Args>
void ReleaseMalloc(T x, Args... args) {
	free(x); ReleaseMalloc(args...);
}
