#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"
#include <vector>

// 如果你不需要使用 matrix 类，请将 IGNORE_MATRIX 改为 0
// #define IGNORE_MATRIX 0
#define IGNORE_MATRIX 0

#if IGNORE_MATRIX

class matrix {
private:
    int m, n;
    fraction **data;
public:
    matrix() { m = n = 0; data = nullptr; }
    matrix(int m_, int n_);
    matrix(const matrix &obj);
    matrix(matrix &&obj) noexcept;
    ~matrix();
    matrix &operator=(const matrix &obj);
    fraction &operator()(int i, int j);
    friend matrix operator*(const matrix &lhs, const matrix &rhs);
    matrix transposition();
    fraction determination();
};

#endif

class resistive_network {
private:
    int n, m;
    struct Edge { int u, v; fraction r; fraction g; };
    std::vector<Edge> edges;
    std::vector<std::vector<fraction>> L; // Laplacian

    // Build reduced (n-1)x(n-1) Laplacian by grounding node n
    void build_reduced(std::vector<std::vector<fraction>> &A) const {
        A.assign(n - 1, std::vector<fraction>(n - 1, fraction(0)));
        for (int i = 0; i < n - 1; ++i) {
            for (int j = 0; j < n - 1; ++j) {
                A[i][j] = L[i][j];
            }
        }
    }

    // Gaussian elimination solver with exact fraction arithmetic
    std::vector<fraction> solve(std::vector<std::vector<fraction>> A, std::vector<fraction> b) const {
        int sz = (int)A.size();
        // Augment
        for (int i = 0; i < sz; ++i) A[i].push_back(b[i]);

        int r = 0;
        for (int c = 0; c < sz && r < sz; ++c) {
            int piv = -1;
            for (int i = r; i < sz; ++i) {
                if (!(A[i][c] == fraction(0))) { piv = i; break; }
            }
            if (piv == -1) continue;
            if (piv != r) std::swap(A[piv], A[r]);
            fraction inv_piv = fraction(1) / A[r][c];
            for (int j = c; j <= sz; ++j) A[r][j] = A[r][j] * inv_piv;
            for (int i = 0; i < sz; ++i) if (i != r) {
                fraction factor = A[i][c];
                if (factor == fraction(0)) continue;
                for (int j = c; j <= sz; ++j) {
                    A[i][j] = A[i][j] - factor * A[r][j];
                }
            }
            ++r;
        }

        std::vector<fraction> x(sz, fraction(0));
        for (int i = 0; i < sz; ++i) x[i] = A[i][sz];
        return x;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        n = interface_size_;
        m = connection_size_;
        edges.resize(m);
        L.assign(n, std::vector<fraction>(n, fraction(0)));
        for (int k = 0; k < m; ++k) {
            int a = from[k];
            int b = to[k];
            fraction r = resistance[k];
            fraction g = fraction(1) / r;
            edges[k] = {a, b, r, g};
            int i = a - 1, j = b - 1;
            L[i][i] = L[i][i] + g;
            L[j][j] = L[j][j] + g;
            L[i][j] = L[i][j] - g;
            L[j][i] = L[j][i] - g;
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int s = interface_id1;
        int t = interface_id2;
        // Build current vector I with +1 at s, -1 at t
        std::vector<fraction> I(n, fraction(0));
        I[s - 1] = I[s - 1] + fraction(1);
        I[t - 1] = I[t - 1] - fraction(1);

        // Reduced system
        std::vector<std::vector<fraction>> A;
        build_reduced(A);
        std::vector<fraction> b(n - 1, fraction(0));
        for (int i = 0; i < n - 1; ++i) b[i] = I[i];
        std::vector<fraction> U = solve(A, b); // potentials for nodes 1..n-1, u_n = 0

        fraction us = (s == n) ? fraction(0) : U[s - 1];
        fraction ut = (t == n) ? fraction(0) : U[t - 1];
        return us - ut; // voltage difference equals effective resistance for 1A
    }

    fraction get_voltage(int id, fraction current[]) {
        // Build reduced system L' * U' = I'
        std::vector<std::vector<fraction>> A;
        build_reduced(A);
        std::vector<fraction> b(n - 1, fraction(0));
        for (int i = 0; i < n - 1; ++i) b[i] = current[i];
        std::vector<fraction> U = solve(A, b);
        return U[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (const auto &e : edges) {
            fraction du = voltage[e.u - 1] - voltage[e.v - 1];
            fraction term = e.g * (du * du);
            total = total + term;
        }
        return total;
    }
};


#endif //SRC_HPP

