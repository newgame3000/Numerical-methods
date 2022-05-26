#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <functional>
#include "polynomial.hpp"


using namespace std;

function<double(double)> f{
    [](double x) { return 1 /((2 * x + 7) * (3 * x + 4)); },
};


double SimpsonMethod(double a, double b, double h) {
    double res = 0;

    int count = floor((b - a) / h);

    double x = a;
    for (int i = 0; i < count; ++i) {
        vector<double> xx {
            x,
            x + h / 2,
            x + h
        };

        vector<double> y {
            f(x),
            f(x + h / 2),
            f(x + h)
        };

        polynomial<double> p;
        p.LagrangePolynomial(xx, y);

        res += p.Integrate(x, x + h);
        x += h;
        //res.Print();
    }

    return res;
}

double TrapezoidalMethod(double a, double b, double h) {
    double res = 0;
    int count = floor((b - a) / h);
    double x = a;

    for (int i = 0; i < count; ++i) {
        res += (f(x) + f(x + h)) * h;
        x += h;
    }
    return res  / 2;
}

double RectangleMethod(double a, double b, double h) {
    double res = 0;
    int count = floor((b - a) / h);
    double x = a;
    for (int i = 0; i < count; ++i) {
        res += f(x) * h;
        x += h;
    }

    return res;
}

double RungeRombergMethod(double I, double r, int m, double a, double b, double h) {
    if (m == 0) {
        double Ir = RectangleMethod(a, b, h * r);
        return abs(abs(Ir - I) / (pow(r, 1) - 1));
    }

    if (m == 1) {
        double Ir = TrapezoidalMethod(a, b, h * r);
        return abs(abs(Ir - I) / (pow(r, 2) - 1));
    }


    if (m == 2) {
        double Ir = SimpsonMethod(a, b, h * r);
        return abs(abs(Ir - I) / (pow(r, 4) - 1));
    }

    return 0;
}


int main() {
 
    double a, b, h;
    cin >> a >> b >> h;

    double res = RectangleMethod(a, b, h);
    cout << "Метод прямоугольников: " << res << endl;
    cout << "Погрешность: " << RungeRombergMethod(res, 0.0001, 0, a ,b, h) << endl;

    res = TrapezoidalMethod(a, b, h);
    cout << "Метод трапеций: " << res << endl;
    cout << "Погрешность: " << RungeRombergMethod(res, 0.0001, 1, a ,b, h) << endl;

    res = SimpsonMethod(a, b, h);
    cout << "Метод Симпсона: " << res << endl;
    cout << "Погрешность: " << RungeRombergMethod(res, 0.0001, 2, a ,b, h) << endl;

}
