#include <iostream>
#include <math.h>

using namespace std;

template <class T>
int Sign(T a){
    if (a > 0) {
        return 1;
    } else {
        if (a < 0) {
            return -1;
        } else {
            return 0;
        }
    }
}

double f(double x) {
    return cos(x) + 0.25 * x - 0.5;
}

double f_1(double x) {
    return -sin(x) + 0.25 ;
}

double Simple_iterations(double l, double r, int & count, double eps) {
    double x;
    double oldx = l;
    double lambda = Sign(f_1(oldx)) / (max(abs(f_1(l)), abs(f_1(r))));
    double q = max(abs(1 - lambda * f_1(l)), abs(1 - lambda * f_1(r)));
    q = q / (1 - q);
    double epsk;

    do {
        count += 1;
        x = oldx - lambda * f(oldx);
        epsk = q * abs(oldx - x);
        oldx = x;
    } while(epsk > eps);
    return oldx;
}

double Newtons_method(double x, int & count, double eps) {
    double oldx = x;
    double epsk;
    do {
        count += 1;
        x = oldx - f(oldx) / f_1(oldx);

        epsk = abs(oldx - x);

        oldx = x;

    } while(epsk > eps);
    return oldx;
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
    double eps;
    cin >> eps;
    double l = -1;
    double r = 0;

    int count = 0;

    cout << fixed;
    int k = DegreeEps(eps);
    cout.precision(k + 1);

    double x  = Simple_iterations(l, r, count, eps);
    cout << "Метод простых итераций\nРешение\n";
    cout << x << endl << "Количество итераций\n" << count << endl;

    count = 0;

    x  = Newtons_method(l, count, eps);
    cout << "Метод Ньютона\nРешение\n";
    cout << x << endl << "Количество итераций\n" << count << endl;

}