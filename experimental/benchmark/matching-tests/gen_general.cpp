/*
gen_general.cpp

Generates a random, undirected graph and evaluates matching to test LAGraph_MaximalMatching

Usage:
Identical to gen_bipartite

Details:
The maximum matching technique used here is an implementation of a lesser known algorithm called the Blossom algorithm.
Read more here: https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/paths-trees-and-flowers/08B492B72322C4130AE800C0610E0E21
The implementation is provided here: https://codeforces.com/blog/entry/92339

The naive technique first sorts all edges by weight for weighted graphs (using the prefer_light option, unweighted graphs have edge weight 1). The sort breaks ties between edges
of the same weight by edge degree (edges with a lower degree are favored). The method then traverses the edge list and attempts to include each edge in order in the matching.

In both cases, the generated graph is printed in MatrixMarket format to data.mtx. This is so the random graph can be evaluated using
GraphBLAS. In addition, the evaluated matching value is printed to stdout (this can be piped out to any file).
*/

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

// #define dbg

char msg[1024];

#define OK(method)                                                  \
{                                                                   \
    int info = method ;                                             \
    if (!(info == GrB_SUCCESS || info != GrB_NO_VALUE))             \
    {                                                               \
        printf ("error! line %d info %d\n", __LINE__, info) ;       \
        printf ("msg is %s\n", msg) ;                               \
        abort ( ) ;                                                 \
    }                                                               \
}

using namespace std;
using ll = long long;

const string HEADER = "%%MatrixMarket matrix coordinate pattern symmetric\n%%GraphBLAS type bool";
const string WEIGHTED_HEADER = "%%MatrixMarket matrix coordinate integer symmetric\n%%GraphBLAS type uint32_t";

bool prefer_light = 0;
bool weighted = 0;
vector<int> deg;

bool cmp_basic_heavy(vector<ll> &a, vector<ll> &b){
    return (a[2] > b[2]);
}

bool cmp_basic_light(vector<ll> &a, vector<ll> &b){
    return (a[2] < b[2]);
}
// a: {u, v, weight}
// u <-> v, weight
bool cmp_with_degree_heavy(vector<ll> &a, vector<ll> &b){
    if(a[2] == b[2]){
        int sum_deg_a = max(deg[a[0]], deg[a[1]]);
        int sum_deg_b = max(deg[b[0]], deg[b[1]]);
        return (sum_deg_a < sum_deg_b);
    }
    return cmp_basic_heavy(a, b);
}

bool cmp_with_degree_light(vector<ll> &a, vector<ll> &b){
    if(a[2] == b[2]){
        int sum_deg_a = max(deg[a[0]], deg[a[1]]);
        int sum_deg_b = max(deg[b[0]], deg[b[1]]);
        return (sum_deg_a < sum_deg_b);
    }
    return cmp_basic_light(a, b);
}

/*
    Credits: Riley Borgard
    https://codeforces.com/blog/entry/92339
*/
struct BlossomSolver {
    int n, m;
    vector<int> mate;
    vector<vector<int>> b;
    vector<int> p, d, bl;
    vector<vector<int>> g;
    BlossomSolver(int n) : n(n) {
        m = n + n / 2;
        mate.assign(n, -1);
        b.resize(m);
        p.resize(m);
        d.resize(m);
        bl.resize(m);
        g.assign(m, vector<int>(m, -1));
    }
    void add_edge(int u, int v) {
        g[u][v] = u;
        g[v][u] = v;
    }
    void match(int u, int v) {
        g[u][v] = g[v][u] = -1;
        mate[u] = v;
        mate[v] = u;
    }
    vector<int> trace(int x) {
        vector<int> vx;
        while(true) {
            while(bl[x] != x) x = bl[x];
            if(!vx.empty() && vx.back() == x) break;
            vx.push_back(x);
            x = p[x];
        }
        return vx;
    }
    void contract(int c, int x, int y, vector<int> &vx, vector<int> &vy) {
        b[c].clear();
        int r = vx.back();
        while(!vx.empty() && !vy.empty() && vx.back() == vy.back()) {
            r = vx.back();
            vx.pop_back();
            vy.pop_back();
        }
        b[c].push_back(r);
        b[c].insert(b[c].end(), vx.rbegin(), vx.rend());
        b[c].insert(b[c].end(), vy.begin(), vy.end());
        for(int i = 0; i <= c; i++) {
            g[c][i] = g[i][c] = -1;
        }
        for(int z : b[c]) {
            bl[z] = c;
            for(int i = 0; i < c; i++) {
                if(g[z][i] != -1) {
                    g[c][i] = z;
                    g[i][c] = g[i][z];
                }
            }
        }
    }
    vector<int> lift(vector<int> &vx) {
        vector<int> A;
        while(vx.size() >= 2) {
            int z = vx.back(); vx.pop_back();
            if(z < n) {
                A.push_back(z);
                continue;
            }
            int w = vx.back();
            int i = (A.size() % 2 == 0 ? find(b[z].begin(), b[z].end(), g[z][w]) - b[z].begin() : 0);
            int j = (A.size() % 2 == 1 ? find(b[z].begin(), b[z].end(), g[z][A.back()]) - b[z].begin() : 0);
            int k = b[z].size();
            int dif = (A.size() % 2 == 0 ? i % 2 == 1 : j % 2 == 0) ? 1 : k - 1;
            while(i != j) {
                vx.push_back(b[z][i]);
                i = (i + dif) % k;
            }
            vx.push_back(b[z][i]);
        }
        return A;
    }
    int solve() {
        for(int ans = 0; ; ans++) {
            fill(d.begin(), d.end(), 0);
            queue<int> Q;
            for(int i = 0; i < m; i++) bl[i] = i;
            for(int i = 0; i < n; i++) {
                if(mate[i] == -1) {
                    Q.push(i);
                    p[i] = i;
                    d[i] = 1;
                }
            }
            int c = n;
            bool aug = false;
            while(!Q.empty() && !aug) {
                int x = Q.front(); Q.pop();
                if(bl[x] != x) continue;
                for(int y = 0; y < c; y++) {
                    if(bl[y] == y && g[x][y] != -1) {
                        if(d[y] == 0) {
                            p[y] = x;
                            d[y] = 2;
                            p[mate[y]] = y;
                            d[mate[y]] = 1;
                            Q.push(mate[y]);
                        }else if(d[y] == 1) {
                            vector<int> vx = trace(x);
                            vector<int> vy = trace(y);
                            if(vx.back() == vy.back()) {
                                contract(c, x, y, vx, vy);
                                Q.push(c);
                                p[c] = p[b[c][0]];
                                d[c] = 1;
                                c++;
                            }else {
                                aug = true;
                                vx.insert(vx.begin(), y);
                                vy.insert(vy.begin(), x);
                                vector<int> A = lift(vx);
                                vector<int> B = lift(vy);
                                A.insert(A.end(), B.rbegin(), B.rend());
                                for(int i = 0; i < (int) A.size(); i += 2) {
                                    match(A[i], A[i + 1]);
                                    if(i + 2 < (int) A.size()) add_edge(A[i + 1], A[i + 2]);
                                }
                            }
                            break;
                        }
                    }
                }
            }
            if(!aug) return ans;
        }
    }
};

int main(int argc, char **argv){
    int num_nodes = atoi(argv[1]);
    double sparse_factor = atof(argv[2]);
    int perf = atoi(argv[3]);

    int naive;
    if(perf){
        naive = 1;
    } else {
        naive = atoi(argv[4]);
    }
    
    if(naive){
        weighted = atoi(argv[5]);
        if(weighted){
            prefer_light = atoi(argv[6]);
        }
    }
    int n = num_nodes;
    if(!naive){
        assert(n <= 1000);
    }

    GrB_Matrix A = NULL;

    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<ll> seed_distr(1, 1e15);
    uint64_t seed = 62 ;// seed_distr(gen);

    GrB_Index *rows, *cols ;
    uint32_t *vals ;

    OK ( LAGraph_Init (msg)) ;

    OK ( LAGraph_Random_Matrix (&A, GrB_UINT32, n, n, (sparse_factor / n), seed, msg)) ;

    GrB_Index nvals ;
    OK ( GrB_Matrix_nvals (&nvals, A)) ;

    OK ( LAGraph_Malloc ((void**)(&rows), nvals, sizeof(GrB_Index), msg)) ;
    OK ( LAGraph_Malloc ((void**)(&cols), nvals, sizeof(GrB_Index), msg)) ;
    OK ( LAGraph_Malloc ((void**)(&vals), nvals, sizeof(uint32_t), msg)) ;

    OK ( GrB_Matrix_extractTuples_UINT32 (rows, cols, vals, &nvals, A)) ;

    deg.resize(n + 1);
    fill(deg.begin(), deg.end(), 0);

    vector<vector<ll>> edges;

    for(GrB_Index i = 0; i < nvals; i++){
        ll u = rows[i] + 1;
        ll v = cols[i] + 1;
        if(u <= v){
            // ignore self-loops, edges above diagonal
            continue;
        }
        ll weight = vals[i];
        if (!weighted){
            weight = 1;
        }
        // cout << "adding edge " << u << " " << v << " " << weight << endl;
        assert(weight >= 0);
        deg[u]++;
        deg[v]++;
        edges.pb(vector<ll> {u, v, weight});
    }

    OK ( LAGraph_Free ((void**)(&rows), msg)) ;
    OK ( LAGraph_Free ((void**)(&cols), msg)) ;
    OK ( LAGraph_Free ((void**)(&vals), msg)) ;
    OK ( GrB_Matrix_free (&A)) ;

    ofstream graph_out("data.mtx");

    if(naive){
        auto start = chrono::high_resolution_clock::now();
        // unordered_set<int> touched;
        
        bool touched[n + 1];
        memset(touched, 0, sizeof(touched));

        uint64_t tot_weight = 0;
        int chosen = 0;

        // if weighted, will sort by both edge weight and degree, if unweighted by degree
        if(prefer_light){
            sort(edges.begin(), edges.end(), cmp_with_degree_light);
        } else {
            sort(edges.begin(), edges.end(), cmp_with_degree_heavy);
        }
        for(auto e : edges){
            int u = e[0]; int v = e[1];
            /*
            if(touched.count(u) || touched.count(v)){
                continue;
            }
            touched.insert(u); touched.insert(v);
            */
            
            if (touched[u] || touched[v]){
                continue;
            }
            touched[u] = touched[v] = 1;
            
            chosen++;
            tot_weight += (weighted ? e[2] : 1);
        }
        if(perf){
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> diff = end - start;
            printf("%.10f\n", diff.count());
        } else {
            cout << tot_weight << endl;
        }
    } else {
        // use blossom algorithm
        BlossomSolver blossom(n);
        for(auto e : edges){
            blossom.add_edge(e[0] - 1, e[1] - 1);
        }
        cout << blossom.solve() << endl;
    }
    if(weighted){
        graph_out << WEIGHTED_HEADER << endl;
    } else {
        graph_out << HEADER << endl;
    }
    graph_out << n << " " << n << " " << edges.size() << endl;
    for(auto elem : edges){
        int u = elem[0];
        int v = elem[1];
        ll w = elem[2];
        if (weighted){
            graph_out << u << " " << v << " " << w << endl;
        } else {
            graph_out << u << " " << v << endl;
        }
    }
}
