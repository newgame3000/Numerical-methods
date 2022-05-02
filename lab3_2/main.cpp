#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include "polynomial.hpp"

using namespace std;

vector<double> Three_diag(int n, vector<double> m, vector<double> d) {
    vector<double> x(n);
    vector<double> q(n);
    vector<double> p(n);

    p[0] = -m[1] / m[0];
    q[0] = d[0] / m[0];
    
    int k = 1;
    for (uint i = 2; i < m.size() - 2; i += 3) { //m[i] - a, m[i + 1] - b, m[i + 2] - c
        p[k] = - m[i + 2] / (m[i + 1] + m[i] * p[k - 1]);
        q[k] = (d[k] - m[i] * q[k - 1])/(m[i + 1] + m[i] * p[k - 1]);
        k += 1;
    }

    p[k] = 0;
    q[k] = (d[k] - m[m.size() - 2] * q[k - 1]) / (m[m.size() - 1] + m[m.size() - 2] * p[k - 1]);

    x[n - 1] = q[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return x;
}

vector<polynomial<double>> Spline(vector<double> x, vector<double> y, int n) {
    if (x.size() != y.size()) {
        return vector<polynomial<double>>();
    } 

    vector<double> h(n - 1);
    for (int i = 1; i < n; ++i) {
        h[i - 1] = x[i] - x[i - 1];
    }

    //n - 2 неизв, так как c_0 = 0
    vector<double> m((n - 4) * 3 + 4);
    vector<double> d2(n - 2);

    m[0] = 2 * (h[0] + h[1]);
    m[1] = h[1];
    d2[0] = 3 * ((y[2] - y[1]) / h[1] - (y[1] - y[0]) / h[0]);
    int j = 2;


    for (int i = 2; i < n - 2; ++i) {
        m[j] = h[i - 2];
        //cout << m[j] << endl;
        m[j + 1] = 2 * (h[i - 2] + h[i - 1]);
        m[j + 2] = h[i - 1];
        d2[i - 1] = 3 * ((y[i + 1] - y[i]) / h[i - 1] - (y[i] - y[i - 1]) / h[i - 2]);
        j += 3;
    }

    m[m.size() - 2] = h[n - 3];
    m[m.size() - 1] = 2 * (h[n - 3] + h[n - 2]);

    d2[n - 3] = 3 * ((y[n - 1] - y[n - 2]) / h[n - 2] - (y[n - 2] - y[n - 3]) / h[n - 3]);

    vector<double> c2 = Three_diag(n - 1, m, d2);

    vector<double> a(n - 1);
    vector<double> b(n - 1);
    vector<double> c(n - 1);
    vector<double> d(n - 1);

    vector<polynomial<double>> spline(n - 1, polynomial<double>(4));

    c[0] = 0;
    for (int i = 0; i < n - 2; ++i) {

        c[i + 1] = c2[i];

        a[i] = y[i];
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3.;
        d[i] = (c[i + 1] - c[i]) / h[i] / 3.;
    }

    a[n - 2] = y[n - 2];
    b[n - 2] = (y[n - 1] - y[n - 2]) / h[n - 2] - 2. / 3. * h[n - 2] * c[n - 2];
    d[n - 2] = - c[n - 2] / (3 * h[n - 2]);

    for (int i = 0; i < n - 1; ++i) {
        polynomial<double> xb(2);
        xb[0] = - x[i];
        xb[1] = 1;
        polynomial<double> xc = xb * xb;
        polynomial<double> xd = xc * xb;
        spline[i] = xd * d[i] +  xc * c[i] +  xb * b[i] + a[i];
    }

    return spline;
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

    vector<polynomial<double>> spline = Spline(x, y, n);

    for (int i = 0; i < n - 1; ++i) {

        cout << "x = [" << x[i] << " " << x[i + 1] << "]\n";

        spline[i].Print();

        cout << endl;
    }
}