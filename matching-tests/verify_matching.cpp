#include <bits/stdc++.h>

#define pb push_back
#define f first
#define s second

using namespace std;

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

pair<int, int> parsePair(string ln){
    pair<int, int> res = {-1, -1};
    string raw = getRaw(ln);
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
    return make_pair(stoi(raw_f), stoi(raw_s));
}

int main(){
    ifstream data_in("grb_result.txt");
    set<int> edges;
    map<int, pair<int, int>> edge_map;
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
            pair<int, int> val = parsePair(s);
            int edge_id = val.s;
            int node_id = val.f;
            if(!edge_map.count(edge_id)){
                edge_map[edge_id] = make_pair(node_id, -1);
            } else {
                if(edge_map[edge_id].s != -1){
                    // bad; more than 2 nodes per edge
                    printf("[ERR] edge (%d) has more than 2 nodes\n", edge_id);
                    return 0;
                }
                edge_map[edge_id].s = node_id;
            }
        }
    }
    map<int, int> freq;
    for(int e : edges){
        if(!edge_map.count(e)){
            printf("[ERR] chosen edge (%d) not in E\n", e);
            return 0;
        }
        auto nodes = edge_map[e];
        freq[nodes.f]++;
        freq[nodes.s]++;
        if(freq[nodes.f] > 1){
            printf("[ERR] node (%d) touched more than once", nodes.f);
            return 0;
        }
        if(freq[nodes.s] > 1){
            printf("[ERR] node (%d) touched more than once", nodes.s);
            return 0;
        }
    }
    printf("Verification passed\n");
}