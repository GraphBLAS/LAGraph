/*
Takes an undirected graph from stdin, a parent mapping, and computes the coarsened graph
Use this to test the results of coarsenings
*/
#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define f first
#define s second

using vvi = vector<vector<int>>;

void transpose(vvi &mat){
    for(int i = 0; i < mat.size(); i++){
        for(int j = 0; j < mat.size(); j++){
            if(i < j){
                swap(mat[i][j], mat[j][i]);
            }
        }
    }
}

vvi mul(vvi a, vvi b){
    assert(a[0].size() == b.size());
    vvi res(a.size(), vector<int>(b[0].size()));
    int n = res.size(); int m = res[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            for(int k = 0; k < a[0].size(); k++){
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return res;
}

int main(){
    // nodes, edges
    int n, m;
    cin >> n >> m;
    vvi adj(n, vector<int>(n));
    vvi s(n, vector<int>(n));
    vvi s_trans;
    for(int i = 0; i < m; i++){
        int u, v;
        cin >> u >> v;
        u--; v--;
        adj[u][v] = 1;
        adj[v][u] = 1;
    }
    for(int i = 0; i < n; i++){
        int par;
        cin >> par;
        par--;
        s[par][i] = 1;
    }
    s_trans = s;
    transpose(s_trans);
    vvi res = mul(s, mul(adj, s_trans));
    for(int i = 0; i < res.size(); i++){
        for(int j = 0; j < res[0].size(); j++){
            cout << res[i][j] << " ";
        }
        cout << endl;
    }
}