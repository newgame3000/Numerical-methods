#include <iostream>
#include "matrix.hpp"
#include <complex>
#include <unistd.h>

using namespace std;


double NormVector(vector<double> x) {
    int sum = 0;
    for (uint i = 0; i < x.size(); ++i) {
        sum += x[i] * x[i];
    }

    return sqrt(sum);
}

vector<double> Calculation_eps1(matrix<double> m, int n) {
    vector<double> sum(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            sum[i] += m[j][i] * m[j][i];
        }
        sum[i] = sqrt(sum[i]);
    }
    return sum;
}


vector<double> Calculation_eps2(matrix<double> m, int n) {
    vector<double> sum(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 2; j < n; ++j) {
            sum[i] += m[j][i] * m[j][i];
        }
        sum[i] = sqrt(sum[i]);
    }
    return sum;
}

pair<complex<double>, complex<double>> Calculation_lambda_complex(matrix<double> m, int j, bool &con, double eps) {
    double re;
    double im;

    double d = (m[j][j] + m[j + 1][j + 1]) * (m[j][j] + m[j + 1][j + 1]) - 4 * (m[j][j] * m[j + 1][j + 1] - m[j][j + 1] * m[j + 1][j]);


    re = (m[j][j] + m[j + 1][j + 1]) / 2;
    im = sqrt(abs(d)) / 2;

    pair<complex<double>, complex<double>> res;

    if ((d > eps) || (im < eps)) {
        con = true;
    }

    res.first = complex<double>(re, im);
    res.second = complex<double>(re, -im);

    return res;

}

void Calculation_lambda(matrix<double> m, vector<complex<double>> &lambda, int n, double eps, bool &con) {
    for (int i = 0; i < n; ++i) {
        if ((i == n - 1) || (abs(m[i + 1][i]) < eps)) {
            lambda[i] = m[i][i];
        } else {            

            pair<complex<double>, complex<double>> k = Calculation_lambda_complex(m, i, con, eps);
            
            if (!con) {
                lambda[i] = k.first;
                i += 1;
                lambda[i] = k.second;
            }
        }
    }
}

void QR_algorithm(matrix<double> m, int n, double eps, vector<complex<double>> &lambda, int &iter) {

    matrix<double> q(n, n);
    matrix<double> r(n, n);

    vector<complex<double>> oldlambda(n, 1e18);

    bool flag12 = true;
    bool flag3 = true;

    while(flag12 || flag3) {

        iter += 1;

        m.QR(q, r);
        m = r * q;


        if ((NormVector(Calculation_eps1(m, n)) < eps) || (NormVector(Calculation_eps2(m, n)) < eps)) {
            flag12 = false;
        }

        bool con = false;

        Calculation_lambda(m, lambda, n, eps, con);

        if (con) {
            continue;
        }

        vector<double> eps3(n);

        for (int i = 0; i < n; ++i) {
            eps3[i] = abs(abs(lambda[i]) - abs(oldlambda[i]));
        }

        if (NormVector(eps3) < eps) {
            flag3 = false;
        }

        oldlambda = lambda;
    }
}

int main(){
    int n;
    cin >> n;

    matrix<double> m(n, n);
    m.Read();

    double eps;
    cin >> eps;

    vector<complex<double>> lambda(n);

    int iter = 0;

    QR_algorithm(m, n, eps, lambda, iter);

    cout << "Собственные значения\n";
    for (int i = 0; i < n; ++i) {
        cout << lambda[i] << " ";
    }
    cout << endl;

    cout << "Количество итераций\n";
    cout << iter << endl;
}