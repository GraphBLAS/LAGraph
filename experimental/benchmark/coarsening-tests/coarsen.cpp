extern "C" {
   #include "LAGraph.h"
   #include "LAGraphX.h"
   #include "GraphBLAS.h"
}

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <chrono>
#include <algorithm>
#include <random>
#include <cassert>
#include <cstring>

#define pb push_back
#define f first
#define s second

using namespace std;
using ll = long long;

pair<int, int> parse_line(string s){
    string ff, ss;
    ff = ss = "";
    bool add_to_f, add_to_s;
    add_to_f = add_to_s = 0;
    for(int i = 0; i < s.size(); i++){
        if(s[i] == '('){
            add_to_f = 1;
        } else if (s[i] == ')'){
            add_to_f = 0;
            add_to_s = 1;
        }
        if(add_to_f){
            ff += s[i];
        }
        if(add_to_s){
            ss += s[i];
        }
    }
    reverse(ss.begin(), ss.end());
    while(ss.back() == ' '){
        ss.pop_back();
    }
    reverse(ss.begin(), ss.end());
    return {stoi(f), stoi(s)};
}

int main(){
    ifstream parent_in("parent.mtx");
    
}