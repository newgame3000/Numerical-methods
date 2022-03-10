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

void SimpleIterations (matrix<double> alpha, vector<double> beta, vector<double> &b, double eps, double epsk, double norm) {
    while(epsk > eps) {
        vector<double> b_0 = b;
        b = beta + alpha * b_0;
        epsk = (norm / (1 - norm)) * NormVector(b - b_0);
    }
}


void Seidel (matrix<double> alpha, vector<double> beta, vector<double> &b, double eps, double epsk, double norm, int n) {
    while(epsk > eps) {
        vector<double> b_0 = b;
        for (int i = 0; i < n; ++i) {
            b[i] = 0;
            for (int j = 0; j < n; ++j) {
                b[i] += alpha[i][j] * b[j];
            }
            b[i] += beta[i];
        }
        epsk = (norm / (1 - norm)) * NormVector(b - b_0);
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

    cout << "Матрица α\n";
    jac.Print();
    cout << "Вектор β\n";
    for (int i = 0; i < n; ++i) {
        cout << beta[i] << " ";
    }
    cout << endl;

    cout << "Норма α\n";
    double norm = jac.Norm();
    cout << norm << endl;

    double eps;
    cin >> eps;

    int k = DegreeEps(eps);
    cout.precision(k + 1);

    b = beta;
    double epsk = eps + 1;
    SimpleIterations(jac, beta, b, eps, epsk, norm);

    cout << "Метод Якоби (простых итераций)\n";

    for (int i = 0; i < n; ++i) {
        cout << b[i] << " ";
    }
    cout << endl;

    b = beta;
    epsk = eps + 1;

    Seidel(jac, beta, b, eps, epsk, norm, n);

    cout << "Метод Зейделя\n";

    for (int i = 0; i < n; ++i) {
        cout << b[i] << " ";
    }
    cout << endl;

    //-8.0 2.0 4.0 5.0 

    //vector<double> debug{-8.0, 2.0, 4.0, 5.0 };


    
}