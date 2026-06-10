#include <bits/stdc++.h>
#define FASTIO ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
using namespace std;
using ll = long long;
using pii = pair<int, int>;
const int N = 1e5;
int fa[N], CT_fa[N], near_CLP[N], up_BN[N], down_BN[N];
vector<int> BN, down_CL[N];
struct Edge {
    int id, to;
};
vector<Edge> edge[N];
int B;
namespace TOP_CLUSTER {
    int cur_CL[N], cur_CL_cnt;
    void add_CL(int u, int v) {
        if (!v)
            v = cur_CL[cur_CL_cnt];
        CT_fa[v] = u; near_CLP[u] = u;
        for (int r = v; r != u; r = fa[r])
            near_CLP[r] = r;
        for (int i = 1; i <= cur_CL_cnt; ++i) {
            int r = cur_CL[i];
            int j;
            up_BN[r] = u, down_BN[r] = v, down_CL[v].push_back(r);
            for (j = r; !near_CLP[j]; j = fa[j]);
            near_CLP[r] = near_CLP[j];
        }
        cur_CL_cnt = 0;
    }
    int ST[N], top, rec_top[N];
    int waiting[N], rec_BN[N];
    inline void CL_partition(int u, int FA) {
        fa[u] = FA;
        rec_top[u] = top;
        for (auto it = edge[u].begin(); it < edge[u].end(); ++it) {
            if (it->to == FA) {
                edge[u].erase(it);
                break;
            }
        }
        waiting[u] = 1;
        int BN_cnt = 0;
        for (auto& v : edge[u]) {
            ST[++top] = v.to;
            CL_partition(v.to, u);
            waiting[u] += waiting[v.to];
            if (rec_BN[v.to]) {
                rec_BN[u] = rec_BN[v.to];
                BN_cnt++;
            }
        }
        if (waiting[u] > B || BN_cnt > 1 || !FA) {
            waiting[u] = 0, rec_BN[u] = u, BN.emplace_back(u);
            for (int i = 0, j = rec_top[u] + 1, cnt = 0, cur_down = 0, v; i <= edge[u].size(); ++i) {
                v = (i == edge[u].size() ? 0 : edge[u][i].to);
                if (cnt + waiting[v] > B || (cur_down && rec_BN[v]) || !v) {
                    for (; (j < rec_top[v] || !v) && j <= top; j++)
                        cur_CL[++cur_CL_cnt] = ST[j];
                    add_CL(u, cur_down);
                    cnt = cur_down = 0;
                }
                cnt += waiting[v];
                if (rec_BN[v])
                    cur_down = rec_BN[N];
            }
        }
        top = rec_top[u];
    }
}
