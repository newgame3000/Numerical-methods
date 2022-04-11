#include <iostream>
#include "matrix.hpp"
#include <vector>
#include <unistd.h>
#include <functional>

using namespace std;

vector<pair<int, int>> sw;
int swaps;
vector<double> x0({0.2, 0.8});

vector<double> a({0.1, 0.5});
vector<double> b({1.5, 2.3});


vector<function<double(double, double)>> f{
    [](double x1, double x2){ return x1 - cos(x2) - 1; },
    [](double x1, double x2){ return x2 - log10(x1 + 1) - 2; }
};

matrix<function<double(double, double)>> J {
    { [](double x1, double x2){ return 1; }, [](double x1, double x2){ return sin(x2); } },
    { [](double x1, double x2){ return - log(10) / (x1 + 1); }, [](double x1, double x2){ return 1; } }
};

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

template <class T>
void SolveLU(vector<T> &x, matrix<T> &l, matrix<T> &u, vector<T> b, int n) {
    for (uint i = 0; i < sw.size(); ++i) {
        swap(b[sw[i].first], b[sw[i].second]);
    }

    vector<T> z(n);
    int k = 1;
    z[0] = b[0];
    for (int i = 1; i < n; ++i) {
        z[i] = b[i];
        for (int j = 0; j < k; ++j) {
            z[i] -= l[i][j] * z[j];
        }
        k += 1;
    }

    k = 1;
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

int DegreeEps(double eps) {
    int k = 0;
    while(eps <= 1) {
        eps *= 10;
        k += 1;
    }
    return k;
}

vector<double> Newtons_method(int & count, double eps) {
    double epsk;
    vector<double> x(x0);
    vector<double> dx(x0);

    do {
        count += 1;

        matrix <double> J_k(J.Columns_count(),J.Columns_count());
        for (uint i = 0; i < J.Columns_count(); ++i) {
            for (uint j = 0; j < J.Columns_count(); ++j) {
                J_k[i][j] = J[i][j](x[0], x[1]);
            }
        }

        vector<double> f_1(f.size());
        for (uint i = 0; i < f.size(); ++i) {
            f_1[i] = -f[i](x[0], x[1]);
        }

        matrix<double> u(J_k);
        matrix<double> l(J.Columns_count(),J.Columns_count());

        sw.clear();
        J_k.LU(l, u, sw, swaps, J.Columns_count());

        SolveLU(dx, l, u, f_1, J.Columns_count());

        x = x + dx;

        epsk = NormVector(dx);
    } while (eps < epsk);
    return x;
}

matrix <double> Inverse(matrix <double> J_k) {
    matrix<double> inv(J_k.Columns_count(), J_k.Columns_count());
    vector<double> x(J_k.Columns_count());

    sw.clear();
    matrix<double> u(J_k);
    matrix<double> l(J_k.Columns_count(),J_k.Columns_count());

    J_k.LU(l, u, sw, swaps, J_k.Columns_count());

    for (uint i = 0; i < J_k.Columns_count(); ++i) {
        b = vector<double>(J_k.Columns_count());
        b[i] = 1;

        SolveLU(x, l, u, b, J.Columns_count());

        inv.Fill_in_the_column(x, i); 
    }
    return inv;
}

matrix<double> CalculateM(matrix<function<double(double, double)>> a, double x, double y) {
    matrix<double> res(a.Columns_count(), a.Columns_count());
    for (size_t i = 0; i < a.Columns_count(); ++i) {
        for (size_t j = 0; j < a.Columns_count(); ++j) {
            res[i][j] = a[i][j](x, y);
        }
    }
    return res;
}

vector<double> CalculateV(vector<function<double(double, double)>> a, double x, double y) {
    vector<double> res(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = f[i](x, y);
    }
    return res;
}

vector<double> Simple_iterations(int & count, double eps) { 
    double epsk;
    vector<double> x({a[0], b[0]});
    vector<double> oldx({a[0], b[0]});

    matrix <double> J_1(J.Columns_count(),J.Columns_count());
    for (uint i = 0; i < J.Columns_count(); ++i) {
        for (uint j = 0; j < J.Columns_count(); ++j) {
            J_1[i][j] = J[i][j](a[0], b[0]);
        }
    }

    matrix <double> E(J.Columns_count(),J.Columns_count());
    for (uint i = 0; i < J.Columns_count(); ++i) {
        E[i][i] = 1;
    }    

    J_1 = Inverse(J_1);


    function<matrix<double>(double, double)> Jphi {
        [&E, &J_1](double x, double y){ return E - J_1 * CalculateM(J, x, y); } 
    };

    double q = Jphi(a[0], b[0]).Norm();
    q = max(q, Jphi(a[1], b[0]).Norm());
    q = max(q, Jphi(a[0], b[1]).Norm());
    q = max(q, Jphi(a[1], b[1]).Norm());

    cout << "q = " << q << endl;

    q = q / (1 - q);

    do {
        count += 1;

        x = oldx - J_1 * CalculateV(f, x[0], x[1]);

        epsk = q * NormVector(x - oldx);
        oldx = x;

    } while (eps < epsk);
    return x;    
}

int main(){
    double eps;
    cin >> eps;

    cout << fixed;
    int k = DegreeEps(eps);
    cout.precision(k + 1);

    vector<double> x(2);
    int count = 0;

    x = Simple_iterations(count, eps);
    cout << "Метод Простых итераций\nРешение\n";
    cout << "x1 = " << x[0] << "\nx2 = " << x[1] << endl << "Количество итераций\n" << count << endl;

    count = 0;

    x = Newtons_method(count, eps);
    cout << "Метод Ньютона\nРешение\n";
    cout << "x1 = " << x[0] << "\nx2 = " << x[1] << endl << "Количество итераций\n" << count << endl;
}