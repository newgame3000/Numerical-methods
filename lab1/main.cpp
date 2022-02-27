#include <iostream>
#include "matrix.hpp"
#include <unistd.h>

using namespace std;

int swaps = 0;
//vector<pair<int, int>> sw;

void SwapVector(vector<double> &b, int i, int j) {
    int k = b[i];
    b[i] = b[j];
    b[j] = k;
}

void Maximum(matrix<double> &u, matrix<double> &l, vector<double> &b, int n, int j) {
    double max = u[j][j];
    int index = j;
    for (int i = j; i < n; ++i) {
        if (u[i][j] > u[j][j] && u[i][j] > max) {
            index = i;
            max = u[i][j];
        }
    }

    if (index != j) {
       // sw.push_back(make_pair(j, index));
        swaps += 1;
        SwapVector(b, j, index);
        l.Swap_rows(j, index);
        l.Swap_columns(j, index);
        u.Swap_rows(j, index);
    }
}

void Solution_Lz(vector<double> &z, matrix<double> l, vector<double> b, int n) {
    int k = 1;
    z[0] = b[0];
    for (int i = 1; i < n; ++i) {
        z[i] = b[i];
        for (int j = 0; j < k; ++j) {
            z[i] -= l[i][j] * z[j];
        }
        k += 1;
    }
}

void LU(matrix<double> &l, matrix<double> &m, int n, vector<double> &b) {
    l[0][0] = 1;
    for (int j = 0; j < n; ++j) {
        l[j][j] = 1;
    }
    int k = 1;
    for (int j = 0; j < n; ++j) {
        Maximum(m, l, b, n, j);
        for (int i = k; i < n; ++i) {
            l[i][j] = m[i][j] / m[j][j];
            for (int t = 0; t < n; ++t) {
                m[i][t] -= m[j][t] * l[i][j];
            }
        }
        k += 1;
    }
}

void Solution_Ux(vector<double> &x, matrix<double> u, vector<double> z, int n) {
    int k = 1;
    x[n - 1] = z[n - 1] / u[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = z [i];
        for (int j = n - 1; j >= n - k; --j) {
            x[i] -= u[i][j] * x[j];
        }
        x[i] /= u[i][i];
        k += 1;
    }
}

int main(){
    int n;
    cin >> n;

    matrix<double> u(n, n);
    u.Read();

    matrix<double> debug(u);

    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }

    matrix<double> l(n, n);

    LU(l, u, n, b);
    cout << "Матрица L\n";
    l.Print();
    cout << "Матрица U\n";
    u.Print();

    //Решаем Lz = b
    vector<double> z(n);
    Solution_Lz(z, l, b, n);
    cout.precision(1);
    cout << "Решение Lz = b\n";
    for (int i = 0; i < n; ++i) {
        cout << z[i] << " ";
    }
    cout << endl;

    //Решаем Ux = z
    vector<double> x(n);
    Solution_Ux(x, u, z, n);
    cout << "Решение системы (Ux = z)\n";
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " ";
    }
    cout << endl;

    //Определитель(det(L) = 1, det(U) = a * b * c *...)
    double det = 1;
    for (int i = 0; i < n; ++i) {
        det *= u[i][i];
    }

    if (swaps % 2 != 0) { //Так как переставляли строки местами
        det *= -1;
    }

    cout << "Определитель\n";
    cout << det << endl;

    //Обратная матрица
    matrix<double> m(n, n);
    b = vector<double>(n);
    b[0] = 1;
    Solution_Lz(z, l, b, n);
    Solution_Ux(x, u, z, n);
    m.Fill_in_the_column(x, 0);

    for (int i = 1; i < n; ++i) {
        b[i - 1] = 0;
        b[i] = 1;

        for (int j = 0; j < n; ++j) {
            cout << b[j] << " ";
        }
        cout << endl;

        Solution_Lz(z, l, b, n);
        Solution_Ux(x, u, z, n);
        m.Fill_in_the_column(x, i);
    }

    // for (uint i = 0; i < sw.size(); ++i) {
    //     m.Swap_columns(sw[i].first, sw[i].second);
    //     // m.Print();
    //     // sleep(1);
    // }

    m.Print();

    matrix<double> inv(debug * m);
    inv.Print();

}