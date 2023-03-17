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

enum Semiring {
    ANY_ONE = 0,
    PLUS_TIMES = 1,
    OR_AND = 2
};


void transpose(vvi &mat){
    for(int i = 0; i < mat.size(); i++){
        for(int j = 0; j < mat.size(); j++){
            if(i < j){
                swap(mat[i][j], mat[j][i]);
            }
        }
    }
}

vvi mul(vvi a, vvi b, Semiring s){
    assert(a[0].size() == b.size());
    vvi res(a.size(), vector<int>(b[0].size()));
    int n = res.size(); int m = res[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            for(int k = 0; k < a[0].size(); k++){
                if(s == PLUS_TIMES){
                    res[i][j] += a[i][k] * b[k][j];
                } else if (s == ANY_ONE || s == OR_AND){
                    res[i][j] = res[i][j] || (a[i][k] && b[k][j]);
                }
            }
        }
    }
    return res;
}

void pr(vvi &mat, char *name = "N/A"){
    int n = mat.size();
    int m = mat[0].size();
    printf("printing %d x %d matrix, name = %s:\n", n, m, name);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
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
    pr(s, "s");
    s_trans = s;
    transpose(s_trans);
    /*
    think of this as follows:
    mapping from orig nodes to new nodes
    check, does some orig node touch a new node? depending on semiring, sums all occurrences
        (summing uses PLUS_SECOND, or PLUS_TIMES). S matrix should be structural only.
        (if not summing, use ANY_ONE).
    
    */
    vvi a = mul(s, adj, PLUS_TIMES);
    pr(a, "s times adj");
    /*
    think of this as follows:
    we have our mapping from orig nodes to new nodes
    s_trans is this: each column says which orig nodes point to the new node of this col
    so, by doing a * s_trans, we are saying: 
    */
    vvi b = mul(a, s_trans, PLUS_TIMES);
    pr(b, "final");
}