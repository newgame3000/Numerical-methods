#include <iostream>
#include "matrix.hpp"
#include <unistd.h>

using namespace std;

int count = 0;

void MaxNotDiagonal(matrix<double> &m, int &maxi, int &maxj) {
    double res = m[0][1];
    maxi = 0;
    maxj = 1;
    for (uint i = 0; i < m.Rows_count(); ++i) {
        for (uint j = i + 1; j < m.Columns_count(); ++j) {
            if (abs(m[i][j]) > abs(res)) {
                res = m[i][j];
                maxi = i;
                maxj = j;
            }
        }
    }
}

void Iterations(matrix<double> &m, matrix<double> &u, double eps, int n) {
    double epsk = eps + 1;
    while (epsk >= eps) {
        count += 1;
        int maxi;
        int maxj;

        MaxNotDiagonal(m, maxi, maxj);

        double phi = 0.5 * atan(double(2 * m[maxi][maxj] / (m[maxi][maxi] - m[maxj][maxj])));

        matrix<double> uk(n, n);

        for (int i = 0; i < n; ++i) {
            uk[i][i] = 1;
        }

        uk[maxi][maxi] = cos(phi);
        uk[maxi][maxj] = -sin(phi);
        uk[maxj][maxi] = sin(phi);
        uk[maxj][maxj] = cos(phi);

        u = u * uk;
        m = uk.Transpose() * m * uk;

        epsk = m.SqrtDiagonal();
    }
}

int DegreeEps(double eps) {
    int k = 0;
    while(eps < 1) {
        eps *= 10;
        k += 1;
    }
    return k;
}

int main(){
    int n;
    cin >> n;

	matrix<double> m(n, n);
    m.Read();
    matrix<double> u(n, n);
    for (int i = 0; i < n; ++i) {
        u[i][i] = 1;
    }
    double eps;
    cin >> eps;

    Iterations(m, u, eps, n);

    int k = DegreeEps(eps);
    cout << fixed << k << endl;
    cout.precision(k + 1);

    cout << "Собственные значения\n";

    for (int i = 0; i < n; ++i) {
        cout << m[i][i] << endl;
    }

    cout << "Собственные векторы\n";

    u.Print();


    cout << "Количество итераций: " << count << endl;
}