#include <iostream>
#include "matrix.hpp"

using namespace std;

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

void LU(matrix<double> &l, matrix<double> &m, int n) {
    l[0][0] = 1;
    int k = 1;
    for (int j = 0; j < n; ++j) {
        l[j][j] = 1;
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

    matrix<double> m(n, n);
    m.Read();

    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }

    matrix<double> l(n, n);

    LU(l, m, n);
    cout << "Матрица L\n";
    l.Print();
    cout << "Матрица U\n";
    m.Print();


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
    Solution_Ux(x, m, z, n);
    cout << "Решение системы (Ux = z)\n";
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " ";
    }
    cout << endl;
}