#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <functional>
#include "polynomial.hpp"

using namespace std;

function<double(double)> f = [](double x){ return log(x); };

int main() {
    int n;
    cin >> n;

    vector<double> x(n);
    vector<double> y(n);

    for (int i = 0; i < n; ++i) {
        cin >> x[i];
        y[i] = f(x[i]);
    }

    double testx;
    cin >> testx;

    polynomial<double> lag(n);
    lag.LagrangePolynomial(x, y);
    cout << "Упрощённый интерполяционный многочлен Лагранжа" << endl;
    lag.Print();

    polynomial<double> newt(n);
    newt.NewtonPolynomial(x, y);
    cout << "Упрощённый интерполяционный многочлен Ньютона" << endl;
    newt.Print();


    cout << "Погрешность в точке " << testx << ":\n";
    double testin = lag(testx);
    testx = f(testx);

    cout << abs(testin - testx) << endl;
}