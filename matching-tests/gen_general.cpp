#include <bits/stdc++.h>

#define pb push_back
#define f first
#define s second

// #define dbg

using namespace std;

const int INF = 1e9;
int NUM_NODES;
int SPARSE_FACTOR;

int main(int argc, char **argv){
    NUM_NODES = atoi(argv[0]);
    SPARSE_FACTOR = atoi(argv[1]);

    int n = NUM_NODES;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distr(1, 2 * SPARSE_FACTOR);

    ofstream graph_out("data.mtx");

    vector<vector<int>> adj(n);

    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            int choose = gen(distr) % SPARSE_FACTOR;
            if(choose == 0){
                // edge chosen
                adj[i].pb(j);
                adj[j].pb(i);
            }
        }
    }
}