#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <functional>


using namespace std;

function<double(double)> y_2{
    [](double z) { return z; },
};

function<double(double, double, double)> z_2{
    [](double x, double y, double z) { return (y - (x - 3) * z) / (x * x - 1); },
};

function<double(double)> q{
    [](double x) { return - 1 / (x * x - 1); },
};

function<double(double)> p{
    [](double x) { return (x - 3) / (x * x - 1); },
};

int a = 0, b = 3;
int y_a = -18, y_b = 0;


vector<double> runge(double h, double m) {
    vector<double> u;

    double xk = a;
    double yk = y_a;
    double zk = m;

    u.push_back(yk);

    int count = floor((b - a) / h);

    for (int i = 0; i < count; ++i) {

        double K1 = h * y_2(zk);
        double L1 = h * z_2(xk, yk, zk);

        double K2 = h * y_2(zk + L1 * 0.5);
        double L2 = h * z_2(xk + h * 0.5, yk + K1 * 0.5, zk + 0.5 * L1);

        double K3 = h * y_2(zk + L2 * 0.5);
        double L3 = h * z_2(xk + h * 0.5, yk + K2 * 0.5, zk + 0.5 * L2);

        double K4 = h * y_2(zk + L3);
        double L4 = h * z_2(xk + h, yk + K3, zk + L3);

        double dy = (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
        double dz = (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;

        xk += h;
        yk += dy;
        zk += dz;

        u.push_back(yk);
    }
    return u;
}

vector<double> threeDiag(int n, vector<double> m, vector<double> d) {
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


vector<double> shooting(double h) {

    int m1 = 0.001;

    int m2 = 0.002;

    vector<double> vecphi1 = runge(h, m1);

    vector<double> vecphi2 = runge(h, m2);

    double phi1 = vecphi1[vecphi1.size() - 1] - y_b; 
    double phi2 = vecphi2[vecphi2.size() - 1] - y_b; 

    while (abs(phi2) > 1e-6) {

        double dphi = (phi2 - phi1) / (m2 - m1);

        double dm = - phi2 / dphi;

        m1 = m2;

        m2 += dm;

        phi1 = phi2;

        vecphi2 = runge(h, m2);

        phi2 = vecphi2[vecphi2.size() - 1] - y_b; 
    }

    return vecphi2;
}

//Краевые условия 3 рода
//alpha * y(a) + beta * y'(a) = y0
//delta * y(b) + gamma * y'(b) = y1
//alpha = 1
//delta = 1

vector<double> fdMethod(double h) {

    int n = floor((b - a) / h) + 1;

    vector <double> m((n - 2) * 3 + 4);

    vector <double> d(n);
    double x = a;

    m[0] = h;
    m[1] = 0;

    d[0] = y_a * h;
    
    int j = 2;
    x += h;
    for (int i = 0; i < n - 2; ++i) {
        m[j] = 1 - h * p (x) * 0.5;
        j += 1;
        m[j] = -2 + h * h * q(x);
        j += 1;
        m[j] = 1 + h * p (x) * 0.5;
        j += 1;

        x += h;
    }

    d[n - 1] = y_b * h;

    m[m.size() - 2] = 0;
    m[m.size() - 1] = h;

    vector<double> u = threeDiag(n, m, d);
    return u;

}

double RungeRombergMethodShooting(vector<double> I, double h, double r) {

    vector<double> Ir = shooting(h * 2);

    double res = 0;

    for (uint i = 0, j = 0; i < Ir.size(); ++i, j+= 2) {

        res = max(res, abs(abs(Ir[i] - I[j]) / (pow(r, 2) - 1)));

    }

    return res;
}

double RungeRombergMethodFd(vector<double> I, double h, double r) {

    vector<double> Ir = fdMethod(h * 2);

    double res = 0;

    for (uint i = 0, j = 0; i < Ir.size(); ++i, j+= 2) {

        res = max(res, abs(abs(Ir[i] - I[j]) / (pow(r, 2) - 1)));

    }

    return res;
}


int main() {
    double h;
    cin >> h;

    cout << "Метод cтрельбы" << endl;

    vector<double> u  = shooting(h);

    FILE *shoot;
    if ((shoot = fopen("shoot.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    double x = a;
    string s = "";
    for (uint i = 0; i < u.size(); ++i) {
        cout << x << " " << u[i] << endl;
        s = "";
        s += to_string(x) + "\n";
        fprintf(shoot, s.c_str());
        s = "";
        s += to_string(u[i]) + "\n";
        fprintf(shoot, s.c_str());
        x += h;
    }

    fclose(shoot);

    cout << "Погрешность: " << RungeRombergMethodShooting(u, h, 2) << endl;;

    cout << "Конечно-разностный метод" << endl;

    u  = fdMethod(h);

    FILE *fd;
    if ((fd = fopen("fd.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    x = a;
    s = "";
    for (uint i = 0; i < u.size(); ++i) {
        cout << x << " " << u[i] << endl;
        s = "";
        s += to_string(x) + "\n";
        fprintf(fd, s.c_str());
        s = "";
        s += to_string(u[i]) + "\n";
        fprintf(fd, s.c_str());
        x += h;
    }

    cout << "Погрешность: " << RungeRombergMethodFd(u, h, 2) << endl;;

    fclose(fd);
}
