#include <bits/stdc++.h>
#define FASTIO ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
using namespace std;
using ll = long long;
using pii = pair<int, int>;
const int P = 998244353;
const int g = 3;
const int gi = 332748118;
ll Pow(ll a, ll b) {
    ll res = 1;
    a %= P;
    for (; b; b >>= 1, a = a * a % P)
        if (b & 1) res = res * a % P;
    return res;
}
void NTT(vector<int>& A, int n, int op) {
    vector<int> R(n, 0);
    for (int i = 0; i < n; i++) R[i] = R[i / 2] / 2 + ((i & 1) ? n / 2 : 0);
    for (int i = 0; i < n; i++)
        if (i < R[i])
            swap(A[i], A[R[i]]);
    for (int i = 2; i <= n; i <<= 1) {
        int g1 = Pow(op == 1 ? g : gi, (P - 1) / i);
        for (int j = 0; j < n; j += i) {
            int gk = 1;
            for (int k = j; k < j + i / 2; k++) {
                int x = A[k], y = 1ll * gk * A[k + i / 2] % P;
                A[k] = (x + y) % P;
                A[k + i / 2] = (x - y + P) % P;
                gk = 1ll * gk * g1 % P;
            }
        }
    }
    if (op != 1) {
        int ni = Pow(n, P - 2);
        for (int i = 0; i < n; i++)
            A[i] = 1ll * A[i] * ni % P;
    }
}
template<typename T>
struct polynomial {
    vector<T> coefficent;
    int length;
    polynomial() : length(0) {}
    polynomial(int len) : coefficent(len, 0), length(len) {}
    polynomial(const vector<T>& v) : coefficent(v), length(v.size()) {}
    polynomial(initializer_list<T> l) : coefficent(l), length(l.size()) {}
    T& operator[](int i) { return coefficent[i]; }
    const T& operator[](int i) const { return coefficent[i]; }
    void resize(int len) {
        coefficent.resize(len, 0);
        length = len;
    }
    void normalize() {
        while (length > 1 && coefficent.back() == 0) {
            coefficent.pop_back();
            length--;
        }
    }
};
template<typename T>
polynomial<T> operator+(const polynomial<T>& a, const polynomial<T>& b) {
    polynomial<T> res(max(a.length, b.length));
    for (int i = 0; i < res.length; i++) {
        T valA = i < a.length ? a[i] : 0;
        T valB = i < b.length ? b[i] : 0;
        res[i] = (valA + valB) % P;
    }
    return res;
}
template<typename T>
polynomial<T> operator-(const polynomial<T>& a, const polynomial<T>& b) {
    polynomial<T> res(max(a.length, b.length));
    for (int i = 0; i < res.length; i++) {
        T valA = i < a.length ? a[i] : 0;
        T valB = i < b.length ? b[i] : 0;
        res[i] = (valA - valB + P) % P;
    }
    return res;
}
template<typename T>
polynomial<T> operator*(const polynomial<T>& a, const polynomial<T>& b) {
    if (a.length == 0 || b.length == 0) return polynomial<T>();
    int target = a.length + b.length - 1;
    int limit = 1;
    while (limit < target) limit <<= 1;
    vector<int> A(limit, 0), B(limit, 0);
    for (int i = 0; i < a.length; i++) A[i] = a[i];
    for (int i = 0; i < b.length; i++) B[i] = b[i];
    NTT(A, limit, 1);
    NTT(B, limit, 1);
    for (int i = 0; i < limit; i++) A[i] = 1ll * A[i] * B[i] % P;
    NTT(A, limit, -1);
    A.resize(target);
    return polynomial<T>(A);
}
template<typename T>
polynomial<T> Inv(polynomial<T> a, int deg) {
    polynomial<T> b({ (T)Pow(a[0], P - 2) });
    for (int len = 2; (len / 2) < deg; len <<= 1) {
        int limit = len * 2;
        vector<int> A(limit, 0), B(limit, 0);
        for (int i = 0; i < min(a.length, len); i++) A[i] = a[i];
        for (int i = 0; i < min(b.length, len); i++) B[i] = b[i];
        NTT(A, limit, 1);
        NTT(B, limit, 1);
        for (int i = 0; i < limit; i++) {
            B[i] = 1ll * B[i] * (2ll - 1ll * A[i] * B[i] % P + P) % P;
        }
        NTT(B, limit, -1);
        b.coefficent.assign(B.begin(), B.begin() + len);
        b.length = len;
    }
    b.resize(deg);
    return b;
}
template<typename T>
polynomial<T> Deriv(const polynomial<T>& a) {
    if (a.length <= 1) return polynomial<T>({ 0 });
    polynomial<T> res(a.length - 1);
    for (int i = 1; i < a.length; i++) {
        res[i - 1] = 1ll * a[i] * i % P;
    }
    return res;
}
template<typename T>
polynomial<T> Integ(const polynomial<T>& a) {
    polynomial<T> res(a.length + 1);
    vector<int> inv(a.length + 1, 1);
    for (int i = 2; i <= a.length; i++) {
        inv[i] = 1ll * (P - P / i) * inv[P % i] % P;
    }
    for (int i = 0; i < a.length; i++) {
        res[i + 1] = 1ll * a[i] * inv[i + 1] % P;
    }
    return res;
}
template<typename T>
polynomial<T> Ln(polynomial<T> a, int deg) {
    a.resize(deg);
    polynomial<T> deriv = Deriv(a);
    polynomial<T> inv = Inv(a, deg);
    polynomial<T> res = deriv * inv;
    res.resize(deg - 1);
    res = Integ(res);
    res.resize(deg);
    return res;
}
template<typename T>
polynomial<T> Exp(polynomial<T> a, int deg) {
    polynomial<T> b({ 1 });
    for (int len = 2; (len / 2) < deg; len <<= 1) {
        polynomial<T> lnb = Ln(b, len);
        polynomial<T> sub(len);
        for (int i = 0; i < len; i++) {
            T valA = i < a.length ? a[i] : 0;
            T valLn = i < lnb.length ? lnb[i] : 0;
            sub[i] = (valA - valLn + P) % P;
        }
        sub[0] = (sub[0] + 1) % P;
        b = b * sub;
        b.resize(len);
    }
    b.resize(deg);
    return b;
}
template<typename T>
polynomial<T> operator/(polynomial<T> a, polynomial<T> b) {
    a.normalize(); b.normalize();
    int n = a.length - 1, m = b.length - 1;
    if (n < m) return polynomial<T>({ 0 });
    polynomial<T> ar = a, br = b;
    reverse(ar.coefficent.begin(), ar.coefficent.begin() + n + 1);
    reverse(br.coefficent.begin(), br.coefficent.begin() + m + 1);
    int deg = n - m + 1;
    polynomial<T> inv_br = Inv(br, deg);
    polynomial<T> q = ar * inv_br;
    q.resize(deg);
    reverse(q.coefficent.begin(), q.coefficent.begin() + deg);
    return q;
}
template<typename T>
polynomial<T> operator%(polynomial<T> a, polynomial<T> b) {
    a.normalize(); b.normalize();
    int n = a.length - 1, m = b.length - 1;
    if (n < m) return a;
    polynomial<T> q = a / b;
    polynomial<T> qb = q * b;
    polynomial<T> r = a - qb;
    r.resize(m);
    r.normalize();
    return r;
}
int main() {
    FASTIO;
    return 0;
}
