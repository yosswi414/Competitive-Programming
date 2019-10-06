/*
*   Date    : Oct 6 2019
*   Author  : yosswi414
*   Note    : Verified at : https://judge.yosupo.jp/problem/unionfind
*/

#include<bits/stdc++.h>
using namespace std;
struct UnionFind {
	vector<int> parent;
	UnionFind(int n) {
		parent.resize(n);
		for (int i = 0; i < n; ++i)parent[i] = i;
	}
	int find(int x) {
		return parent[x] == x ? x : (parent[x] = find(parent[x]));
	}
	void unite(int x, int y) {
		if (find(x) != find(y)) parent[find(y)] = find(x);
	}
	bool issame(int x, int y) {
		return find(x) == find(y);
	}
};

int main() {
	int n, q;
	cin >> n >> q;
	UnionFind U(n);
	for (int i = 0; i < q; ++i) {
		int t, u, v;
		cin >> t >> u >> v;
		if (t == 0)U.unite(u, v);
		else cout << U.issame(u, v) << "\n";
	}
}
