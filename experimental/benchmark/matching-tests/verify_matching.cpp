/*
verify_matching.cpp

Usage:
./verify_matching.cpp

Details:
Reads and parses output from the file 'grb_result.txt'. File must contain ONLY the following data in order:
1. Printout of a matching vector (result of LAGraph_MaximalMatching), produced by LAGraph_Vector_Print() with LAGraph_COMPLETE
2. Prinout of E matrix (input of LAGraph_MaximalMatching), produced by LAGraph_Matrix_Print() with LAGraph_COMPLETE

Determines if the matching described in 'grb_result.txt' is a valid matching, and if so returns the matching value (#of edges/sum of edge weights)

Known bugs: If the edge weights are large, LAGraph_Matrix_Print uses scientific notation, which messes up the parsing
*/

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>

#define pb push_back
#define f first
#define s second

using namespace std;
using ll = long long;

string getRaw(string ln){
    string raw = "";
    bool go = 0;
    for(int i = 0; i < ln.size(); i++){
        if(ln[i] == ')'){
            break;
        }
        if(go){
            raw += ln[i];
        }
        if(ln[i] == '('){
            go = 1;
        }
    }
    return raw;
}

int parseSingle(string ln){
    return stoi(getRaw(ln));
}

vector<int> parsePair(string ln){
    vector<int> res(3, -1);
    string raw = getRaw(ln);
    const int FIRST_INDENT = 4;
    const int SECOND_INDENT = 3;
    int at = FIRST_INDENT + 2 + raw.size() + SECOND_INDENT;
    string weight_raw = ln.substr(at, ln.size() - at);
    int w = stoi(weight_raw);
    string raw_f = ""; string raw_s = "";
    int comma_pos = -1;
    for(int i = 0; i < raw.size(); i++){
        if(raw[i] == ','){
            comma_pos = i;
            break;
        }
    }
    for(int i = 0; i < comma_pos; i++){
        raw_f += raw[i];
    }
    // skip the space following comma
    for(int i = comma_pos + 2; i < raw.size(); i++){
        raw_s += raw[i];
    }
    res[0] = stoi(raw_f);
    res[1] = stoi(raw_s);
    res[2] = w;
    return res;
}

int main(){
    ifstream data_in("grb_result.txt");
    set<int> edges;
    map<int, pair<int, int>> edge_map;
    map<int, int> weight_map;
    int which_mat = -1;
    while(!data_in.eof()){
        string s;
        getline(data_in, s);
        if(s.size() == 0){
            continue; // last line
        }
        if(s.substr(0, 3) == "GrB"){
            which_mat++;
            continue;
        }
        if(which_mat == 0){
            int val = parseSingle(s);
            edges.insert(val);
        } else if (which_mat == 1){
            auto val = parsePair(s);
            int edge_id = val[1];
            int node_id = val[0];
            if(weight_map.count(edge_id)){
                if(weight_map[edge_id] != val[2]){
                    // bad; edge has multiple weights
                    cout << -1 << endl;
                    return 0;
                }
            }
            weight_map[edge_id] = val[2];
            if(!edge_map.count(edge_id)){
                edge_map[edge_id] = make_pair(node_id, -1);
            } else {
                if(edge_map[edge_id].s != -1){
                    // bad; more than 2 nodes per edge
                    cout << -1 << endl;
                    return 0;
                }
                edge_map[edge_id].s = node_id;
            }
        }
    }
    map<int, int> freq;
    ll tot_weight = 0;
    for(int e : edges){
        if(!edge_map.count(e)){
            cout << -1 << endl;
            return 0;
        }
        bool not_found = !edge_map.count(e);
        auto nodes = edge_map[e];
        if(not_found || (nodes.f == -1 || nodes.s == -1)){
            cout << -1 << endl;
            return 0;
        }
        freq[nodes.f]++;
        freq[nodes.s]++;
        if(freq[nodes.f] > 1){
            cout << -1 << endl;
            return 0;
        }
        if(freq[nodes.s] > 1){
            cout << -1 << endl;
            return 0;
        }
        tot_weight += weight_map[e];
    }
    cout << tot_weight << endl;
}