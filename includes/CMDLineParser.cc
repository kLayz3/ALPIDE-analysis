#include <cstdio>
#include <regex>
#include <string>
#include <iostream>

using namespace std;

bool IsHelpArg(int argc, char** argv) {
    for(int i(1); i<argc; ++i) {
        if(!strcmp(argv[i], "help") || !strcmp(argv[i], "--help")) return 1;
    }   
    return 0;
}

/* Cmd args have to be: --tag0=identifier0 --tag1=identifier1 ... */
bool ParseCmdLine(const char* line, string& parsed, int argc, char** argv) {
    cmatch m;
    std::regex r("^(--|)([^=]+)[=](.+)$");
    for(int i(1); i<argc; ++i) {
        if(regex_match(argv[i], m, r) && !strcmp(m[2].str().c_str(), line)) {
            parsed = m[3].str(); return 1;
		}
    }   
    return 0;
}




