#include <iostream>
#include "matrix.hpp"

using namespace std;

vector<double> operator+(const vector<double> &l, const vector<double> &r) {
    vector<double> res = l;
    for (uint i = 0; i < l.size(); ++i) {
        res[i] += r[i];
    }
    return res;
}

vector<double> operator-(const vector<double> &l, const vector<double> &r) {
    vector<double> res = l;
    for (uint i = 0; i < l.size(); ++i) {
        res[i] -= r[i];
    }
    return res;
}

double NormVector(vector<double> x) {
    double max = abs(x[0]);
    for (uint i = 0; i < x.size(); ++i) {
        if (max < abs(x[i])) {
            max = abs(x[i]);
        }
    }
    return max;
}

void Jacobi(matrix<double> m, vector<double> &b, matrix<double>& jac, vector<double> &beta, int n) {
    for (int i = 0; i < n; ++i) {
        beta[i] = b[i] / m[i][i];
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                jac[i][j] = -m[i][j] / m[i][i];
            }
        }
    }
}

void SimpleIterations (matrix<double> alpha, vector<double> beta, vector<double> &b, double eps, double epsk, double norm, int &count) {
    while(epsk > eps) {
        count += 1;
        vector<double> b_0 = b;
        b = beta + alpha * b_0;
        if (norm < 1) {
            epsk = (norm / (1 - norm)) * NormVector(b - b_0);
        } else {
            epsk = NormVector(b - b_0);
        }
    }
}


void Seidel (matrix<double> alpha, vector<double> beta, vector<double> &b, double eps, double epsk, double norm, int n, int &count, double cnorm) {
    while(epsk > eps) {
        count += 1;
        vector<double> b_0 = b;
        for (int i = 0; i < n; ++i) {
            b[i] = 0;
            for (int j = 0; j < n; ++j) {
                b[i] += alpha[i][j] * b[j];
            }
            b[i] += beta[i];
        }
        if (norm < 1) {
            epsk = (cnorm / (1 - norm)) * NormVector(b - b_0);
        } else {
            epsk = NormVector(b - b_0);
        }   
    }
}

int DegreeEps(double eps) {
    int k = 0;
    while(eps <= 1) {
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

    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }

    matrix<double> jac(n, n);
    vector<double> beta(n);

    Jacobi(m, b, jac, beta, n);

    cout << "?????????????? ??\n";
    jac.Print();
    cout << "???????????? ??\n";
    for (int i = 0; i < n; ++i) {
        cout << beta[i] << " ";
    }
    cout << endl;

    cout << "?????????? ??\n";
    double norm = jac.Norm();
    cout << norm << endl;

    double eps;
    cin >> eps;

    int k = DegreeEps(eps);
    cout.precision(k + 1);

    b = beta;
    double epsk = eps + 1;

    int count = 0;
    SimpleIterations(jac, beta, b, eps, epsk, norm, count);

    cout << "?????????? ?????????? (?????????????? ????????????????)\n";

    for (int i = 0; i < n; ++i) {
        cout << b[i] << " ";
    }
    cout << endl;

    cout << "???????????????????? ????????????????: " << count << endl;

    b = beta;
    epsk = eps + 1;

    count = 0;

    matrix<double> c = jac;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            c[i][j] = 0;
        }
    }

    double cnorm = c.Norm();

    Seidel(jac, beta, b, eps, epsk, norm, n, count, cnorm);

    cout << "?????????? ??????????????\n";

    for (int i = 0; i < n; ++i) {
        cout << b[i] << " ";
    }
    cout << endl;  

    cout << "???????????????????? ????????????????: " << count << endl;
}