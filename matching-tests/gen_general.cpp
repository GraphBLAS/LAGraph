#include <bits/stdc++.h>

#define pb push_back
#define f first
#define s second

// #define dbg

using namespace std;

const string HEADER = "%%MatrixMarket matrix coordinate integer symmetric\n%%GraphBLAS type uint32_t";

const int INF = 1e9;
bool prefer_light;
map<int, int> deg;

bool cmp_basic(vector<int> &a, vector<int> &b){
    if(prefer_light){
        if(a[2] < b[2]){
            return 1;
        }        
    } else {
        if(a[2] > b[2]){
            return 1;
        }
    }
    return 0;
}

bool cmp_with_degree(vector<int> &a, vector<int> &b){
    if(a[2] == b[2]){
        int sum_deg_a = max(deg[a[0]], deg[a[1]]);
        int sum_deg_b = max(deg[b[0]], deg[b[1]]);
        if(sum_deg_a < sum_deg_b){
            return 1;
        } else {
            return 0;
        }
    }
    return cmp_basic(a, b);
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
    int sparse_factor = atoi(argv[2]);
    bool naive = atoi(argv[3]); // if selected, we just use the naive greedy way
    bool weighted;
    if(naive){
        weighted = atoi(argv[4]);
        if(weighted){
            prefer_light = atoi(argv[5]);
        }
    }
    int n = num_nodes;

    if(!naive){
        assert(n <= 1000);
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distr(1, 2 * sparse_factor);
    uniform_int_distribution<int> weight_distr(1, 1000);

    ofstream graph_out("data.mtx");
    vector<vector<int>> edges;

    for(int i = 1; i <= n; i++){
        for(int j = i + 1; j <= n; j++){
            int choose = distr(gen) % sparse_factor;
            int weight = (naive && weighted) ? weight_distr(gen) : 1;
            if(choose == 0){
                // edge chosen
                deg[i]++;
                deg[j]++;
                edges.pb(vector<int>{i, j, weight});
            }
        }
    }
    if(naive){
        if(!weighted){
            random_shuffle(edges.begin(), edges.end()); // choose a random order
        } else {
            sort(edges.begin(), edges.end(), cmp_basic); // greedily choose by edge weight
        }
        set<int> touched;
        uint64_t tot_weight = 0;
        int chosen = 0;
        /*
        for(auto e : edges){
            int u = e[0]; int v = e[1];
            if(touched.count(u) || touched.count(v)){
                continue;
            }
            touched.insert(u); touched.insert(v);
            chosen++;
            tot_weight += e[2];
        }
        // if weighted, sorts just by edge weight, if unweighted just a random order
        cout << tot_weight << endl;
        touched.clear();
        tot_weight = chosen = 0;
        */
        // smart naive (it beats GraphBLAS ?!?!?)
        // if weighted, will sort by both edge weight and degree, if unweighted by degree
        sort(edges.begin(), edges.end(), cmp_with_degree);
        for(auto e : edges){
            int u = e[0]; int v = e[1];
            if(touched.count(u) || touched.count(v)){
                continue;
            }
            touched.insert(u); touched.insert(v);
            chosen++;
            tot_weight += e[2];
        }
        cout << tot_weight << endl;
    } else {
        // use blossom algorithm
        BlossomSolver blossom(n);
        for(auto e : edges){
            blossom.add_edge(e[0] - 1, e[1] - 1);
        }
        cout << blossom.solve() << endl;
    }
    graph_out << HEADER << endl;
    graph_out << n << " " << n << " " << edges.size() << endl;
    for(auto elem : edges){
        int u = elem[0];
        int v = elem[1];
        int w = elem[2];
        graph_out << u << " " << v << " " << w << endl;
    }

}