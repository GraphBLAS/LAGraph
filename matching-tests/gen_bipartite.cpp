#include <bits/stdc++.h>

#define pb push_back
#define f first
#define s second

// #define dbg

using namespace std;

const int INF = 1e9;
int NUM_NODES;
int SPARSE_FACTOR = 200; // increase this for sparser graphs. Expected #of edges is (n * n) / SPARSE_FACTOR

const string HEADER = "%%MatrixMarket matrix coordinate pattern symmetric\n%%GraphBLAS type bool";

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

int main(int argc, char **argv){
    NUM_NODES = atoi(argv[1]);
    SPARSE_FACTOR = atoi(argv[2]);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distr(1, 2 * SPARSE_FACTOR);

    ofstream graph_out("data.mtx");

    int n = NUM_NODES / 2;
    int m = NUM_NODES / 2;   
    // cout << n << " " << m << endl;
    adj.resize(n + m + 2);
    capacity.resize(n + m + 2);
    for(int i = 0; i < n + m + 2; i++){
        capacity[i].resize(n + m + 2);
        fill(capacity[i].begin(), capacity[i].end(), 0);
    }

    vector<pair<int, int>> edges;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= m; j++){
            int choose = distr(gen) % SPARSE_FACTOR;
            if(choose == 0){
                // add edge
                adj[i].pb(j + n);
                adj[j + n].pb(i);
                capacity[i][j + n] = 1;
                edges.pb({i, j + n});
            }
        }
    }
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
    graph_out << HEADER << endl;
    graph_out << (n + m) << " " << (n + m) << " " << edges.size() << endl;
    for(auto elem : edges){
        int u = elem.f;
        int v = elem.s;
        graph_out << u << " " << v << endl;
    }
    cout << maxflow(0, n + m + 1, n + m + 2) << endl;
}
