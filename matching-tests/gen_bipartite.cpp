#include <bits/stdc++.h>

#define pb push_back
#define f first
#define s second

// #define dbg

using namespace std;
using ll = long long;

const ll INF = 1e18;

const string HEADER = "%%MatrixMarket matrix coordinate integer symmetric\n%%GraphBLAS type uint32_t";

vector<vector<ll>> cost;
int n, m;

/*
Credits: Andrey Lopatin (https://zafar.cc/2017/7/19/hungarian-algorithm/)
*/
void hungarian(ll &C) {
    vector<ll> u (n+1), v (m+1), p (m+1), way (m+1);
    vector<vector<ll>> A = cost;
    for (int i=1; i<=n; ++i) {
        p[0] = i;
        int j0 = 0;
        vector<ll> minv (m+1, INF);
        vector<char> used (m+1, false);
        do {
            used[j0] = true;
            int i0 = p[j0],  j1;
            ll delta = INF;
            for (int j=1; j<=m; ++j)
                if (!used[j]) {
                    ll cur = A[i0][j]-u[i0]-v[j];
                    if (cur < minv[j])
                        minv[j] = cur,  way[j] = j0;
                    if (minv[j] < delta)
                        delta = minv[j],  j1 = j;
                }
            for (int j=0; j<=m; ++j)
                if (used[j])
                    u[p[j]] += delta,  v[j] -= delta;
                else
                    minv[j] -= delta;
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    C = -v[0];
}

// ***
vector<vector<int>> adj;
vector<vector<int>> capacity;

int bfs(int s, int t, vector<int>& parent) {
    fill(parent.begin(), parent.end(), -1);
    parent[s] = -2;
    queue<pair<int, int>> q;
    q.push({s, INF});
    int ans = 0;
    while (!q.empty()) {
        int cur = q.front().first;
        int flow = q.front().second;
        q.pop();

        for (int next : adj[cur]) {
            if (parent[next] == -1 && capacity[cur][next]) {
                parent[next] = cur;
                int new_flow = min(flow, capacity[cur][next]);
                if (next == t)
                    ans = new_flow;
                q.push({next, new_flow});
            }
        }
    }

    return ans;
}

int maxflow(int s, int t, int n) {
    int flow = 0;
    vector<int> parent(n);
    int new_flow;

    while (new_flow = bfs(s, t, parent)) {
        flow += new_flow;
        int cur = t;
        #ifdef dbg
            cout << "chose path: " << endl;
        #endif
        vector<int> verts;
        verts.pb(t);
        while (cur != s) {
            int prev = parent[cur];
            verts.pb(prev);
            capacity[prev][cur] -= new_flow;
            capacity[cur][prev] += new_flow;
            cur = prev;
        }
        reverse(verts.begin(), verts.end());
        #ifdef dbg
            for(int vv : verts){
                cout << vv << " ";
            }   
            cout << endl;
        #endif
    }

    return flow;
}
// ***

int main(int argc, char **argv){
    int num_nodes = atoi(argv[1]);
    int sparse_factor = atoi(argv[2]);
    // bool weighted = atoi(argv[3]);
    bool weighted = 0;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distr(1, 2 * sparse_factor);
    uniform_int_distribution<int> weight_distr(1, 1e9);

    ofstream graph_out("data.mtx");
    assert(num_nodes % 2 == 0);
    assert(num_nodes <= 1e3);
    n = m = num_nodes / 2;
    // cout << n << " " << m << endl;
    // ***
    adj.resize(n + m + 2);
    capacity.resize(n + m + 2);
    for(int i = 0; i < n + m + 2; i++){
        capacity[i].resize(n + m + 2);
        fill(capacity[i].begin(), capacity[i].end(), 0);
    }
    // ***
    cost.resize(n + 1);
    for(int i = 1; i <= n; i++){
        cost[i].resize(m + 1);
    }
    vector<vector<int>> edges;
    // ***
    vector<pair<int, int>> edges2;
    // ***
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= m; j++){
            int choose = distr(gen) % sparse_factor;
            int weight = weighted ? weight_distr(gen) : 1;
            if(choose == 0){
                // add edge
                // one entry in the cost matrix (Hungarian method)
                cost[i][j] = -weight;
                edges.pb(vector<int> {i, j + n, weight});
                // ***
                adj[i].pb(j + n);
                adj[j + n].pb(i);
                capacity[i][j + n] = 1;
                edges2.pb({i, j + n});
                // ***
            } else {
                cost[i][j] = 1;
            }
        }
    }
    cout << edges.size() << endl;
    // ***
    for(int i = 1; i <= n; i++){
        // source
        adj[0].pb(i);
        adj[i].pb(0);
        capacity[0][i] = 1;
    }
    for(int i = n + 1; i <= n + m; i++){
        // sink
        adj[i].pb(n + m + 1);
        adj[n + m + 1].pb(i);
        capacity[i][n + m + 1] = 1;
    }
    // ***

    graph_out << HEADER << endl;
    graph_out << (n + m) << " " << (n + m) << " " << edges.size() << endl;
    for(auto elem : edges){
        int u = elem[0];
        int v = elem[1];
        int w = elem[2];
        graph_out << u << " " << v << " " << w << endl;
    }
    ll C;
    hungarian(C);
    cout << C << endl;
    cout << "mf: " << maxflow(0, n + m + 1, n + m + 2) << endl;
}
