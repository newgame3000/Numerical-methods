#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <functional>
#include "polynomial.hpp"
#include "matrix.hpp"


using namespace std;

vector<pair<int, int>> sw;
int swaps;

vector<function<double(double)>> f{
    [](double) { return 1; },
    [](double x) { return x; },
    [](double x) { return x*x; },
    [](double x) { return x*x*x; }
};

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

polynomial<double> LSM(vector<double> x, vector<double> y, int p) {

    if (x.size() != y.size()) {
        return polynomial<double>();
    } 


    matrix<double> phi(x.size(), p + 1);
    for (uint i = 0; i < x.size(); ++i) {
        for (int j = 0; j < p + 1; ++j) {
            phi[i][j] = f[j](x[i]);
        }
    }

    matrix<double> phit = phi.Transpose();
    matrix<double> g = phit * phi;
    vector<double> z = phit * y;


    matrix<double> u(g);
    matrix<double> l(u.Columns_count(), u.Columns_count());

    sw.clear();
    g.LU(l, u, sw, swaps, g.Columns_count());

    vector<double> dx(p + 1);

    SolveLU(dx, l, u, z, g.Columns_count());

    return(polynomial<double>(dx));
}

double SquareError(vector<double> x, vector<double> y, polynomial<double> p) {
    double res = 0;

    for (uint i = 0; i < x.size(); ++i) {
        res += (y[i] - p(x[i])) * (y[i] - p(x[i]));
    }
    return res;
}

int main() {
    int n;
    cin >> n;

    vector<double> x(n);
    vector<double> y(n);

    for (int i = 0; i < n; ++i) {
        cin >> x[i];
        cin >> y[i];
    }

    cout << "Приближающий многочлен 1-ой степени:\n"; 

    polynomial<double> res = LSM(x, y, 1);
    res.Print();

    cout << "Ошибка\n" << SquareError(x, y, res) << endl;


    cout << "Приближающий многочлен 2-ой степени:\n"; 

    res = LSM(x, y, 2);
    res.Print();

    cout << "Ошибка\n" << SquareError(x, y, res) << endl;

    cout << "Приближающий многочлен 3-ей степени:\n"; 

    res = LSM(x, y, 3);
    res.Print();

    cout << "Ошибка\n" << SquareError(x, y, res) << endl;
}