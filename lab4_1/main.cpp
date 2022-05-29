#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <functional>


using namespace std;

function<double(double)> y_2{
    [](double z) { return z; },
};

function<double(double, double)> z_2{
    [](double x, double y) { return (1 + 2 * tan(x) * tan(x)) * y ; },
};


int a = 0, b = 1;
int y_0 = 1, z_0 = 2;

vector<double> euler(double h) {
    vector<double> u;

    double xk = a;
    double yk = y_0;
    double zk = z_0;

    u.push_back(yk);

    int count = floor((b - a) / h);

    for (int i = 0; i < count; ++i) {

        double dy = h * y_2(zk);
        double dz = h * z_2(xk, yk);

        yk += h * (y_2(zk) + y_2(zk + dy)) / 2;
        zk += h * (z_2(xk, yk) + z_2(xk + h, yk + dz)) / 2;

        xk += h;

        u.push_back(yk);
    }
    return u;
}

vector<double> runge(double h, vector<double> &z) {
    vector<double> u;

    double xk = a;
    double yk = y_0;
    double zk = z_0;

    u.push_back(yk);
    z.push_back(zk);

    int count = floor((b - a) / h);

    for (int i = 0; i < count; ++i) {

        double K1 = h * y_2(zk);
        double L1 = h * z_2(xk, yk);

        double K2 = h * y_2(zk + L1 / 2);
        double L2 = h * z_2(xk + h / 2, yk + K1 / 2);

        double K3 = h * y_2(zk + L2 / 2);
        double L3 = h * z_2(xk + h / 2, yk + K2 / 2);

        double K4 = h * y_2(zk + L3);
        double L4 = h * z_2(xk + h, yk + K3);

        double dy = (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
        double dz = (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;


        xk += h;
        yk += dy;
        zk += dz;

        u.push_back(yk);
        z.push_back(zk);
    }
    return u;
}


vector<vector <double>> adams(double h, vector<vector <double>> first_four) {

    double xk = first_four[first_four.size() - 1][0];
    double yk = first_four[first_four.size() - 1][1];
    double zk = first_four[first_four.size() - 1][2];


    int count = floor((b - a) / h) - 3;

    for (int i = 0; i < count; ++i) {

        double dy = (h / 24.0) * (55.0 * y_2(first_four[first_four.size() - 1][2])
                                    - 59.0 * y_2(first_four[first_four.size() - 2][2])
                                    + 37.0 * y_2(first_four[first_four.size() - 3][2])
                                    - 9.0 * y_2(first_four[first_four.size() - 4][2]));

        double dz = (h / 24.0) * (55.0 * z_2(first_four[first_four.size() - 1][0], first_four[first_four.size() - 1][1])
                                    - 59.0 * z_2(first_four[first_four.size() - 2][0], first_four[first_four.size() - 2][1])
                                    + 37.0 * z_2(first_four[first_four.size() - 3][0], first_four[first_four.size() - 3][1])
                                    - 9.0 * z_2(first_four[first_four.size() - 4][0], first_four[first_four.size() - 4][1]));


        double xk1 = xk + h;
        double yk1 = yk + dy;
        double zk1 = zk + dz;


        first_four.push_back(vector<double>(3));
        first_four[first_four.size() - 1][0] = xk1;
        first_four[first_four.size() - 1][1] = yk1;
        first_four[first_four.size() - 1][2] = zk1;

        dy = (h / 24.0) * (9.0 * y_2(first_four[first_four.size() - 1][2])
                        + 19.0 * y_2(first_four[first_four.size() - 2][2])
                        - 5.0 * y_2(first_four[first_four.size() - 3][2])
                        + 1.0 * y_2(first_four[first_four.size() - 4][2]));


        dz = (h / 24.0) * (9.0 * z_2(first_four[first_four.size() - 1][0], first_four[first_four.size() - 1][1])
                        + 19.0 * z_2(first_four[first_four.size() - 2][0], first_four[first_four.size() - 2][1])
                        - 5.0 * z_2(first_four[first_four.size() - 3][0], first_four[first_four.size() - 3][1])
                        + 1.0 * z_2(first_four[first_four.size() - 4][0], first_four[first_four.size() - 4][1]));
            
        xk += h;
        yk += dy;
        zk += dz;


        first_four.pop_back();
        first_four.push_back(vector<double>(3));

        first_four[first_four.size() - 1][0] = xk;
        first_four[first_four.size() - 1][1] = yk;
        first_four[first_four.size() - 1][2] = zk;
    }

    return first_four;
}


double RungeRombergMethodEuler(vector<double> I, double h, double r) {

    vector<double> Ir = euler(h * 2);

    double res = 0;

    for (uint i = 0, j = 0; i < Ir.size(); ++i, j+= 2) {

        res = max(res, abs(abs(Ir[i] - I[j]) / (pow(r, 2) - 1)));

    }

    return res;
}

double RungeRombergMethodRunge(vector<double> I, double h, double r) {
    vector <double> z;
    vector<double> Ir = runge(h * 2, z);

    double res = 0;

    for (uint i = 0, j = 0; i < Ir.size(); ++i, j+= 2) {

        res = max(res, abs(abs(Ir[i] - I[j]) / (pow(r, 4) - 1)));

    }

    return res;
}

double RungeRombergMethodAdams(vector<vector <double>> I, double h, double r) {

    vector<double> z;

    vector<double> u  = runge(h * 2, z);

    vector<vector <double>> first_four(4, vector<double> (3));

    double x = a;
    for (int i = 0; i < 4; ++i) {
        first_four[i][0] = x;
        first_four[i][1] = u[i];
        first_four[i][2] = z[i];

        x += h * 2;
    }


    vector<vector <double>> Ir = adams(h * 2, first_four);


    double res = 0;

    for (uint i = 0, j = 0; i < Ir.size(); ++i, j+= 2) {

        res = max(res, abs(abs(Ir[i][1] - I[j][1]) / (pow(r, 4) - 1)));

    }

    return res;
}


int main() {
    double h;
    cin >> h;

    cout << "Метод Эйлера" << endl;

    vector<double> u  = euler(h);

    FILE *eul;
    if ((eul = fopen("eul.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    double x = a;
    string s = "";
    for (uint i = 0; i < u.size(); ++i) {
        cout << x << " " << u[i] << endl;
        s = "";
        s += to_string(x) + "\n";
        fprintf(eul, s.c_str());
        s = "";
        s += to_string(u[i]) + "\n";
        fprintf(eul, s.c_str());
        x += h;
    }

    cout << "Погрешность: " << RungeRombergMethodEuler(u, h, 2) << endl;

    fclose(eul);

    vector<double> z;

    cout << "Метод Рунге-Кутта" << endl;

    u  = runge(h, z);

    FILE *runge;
    if ((runge = fopen("runge.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    x = a;
    for (uint i = 0; i < u.size(); ++i) {
        cout << x << " " << u[i] << endl;
        s = "";
        s += to_string(x) + "\n";
        fprintf(runge, s.c_str());
        s = "";
        s += to_string(u[i]) + "\n";
        fprintf(runge, s.c_str());
        x += h;
    }

    cout << "Погрешность: " << RungeRombergMethodRunge(u, h, 2) << endl;

    fclose(runge);

    cout << "Метод Адамса" << endl;

    vector<vector <double>> first_four(4, vector<double> (3));

    x = a;
    for (int i = 0; i < 4; ++i) {
        first_four[i][0] = x;
        first_four[i][1] = u[i];
        first_four[i][2] = z[i];

        x += h;
    }

    first_four = adams(h, first_four);

    FILE *adams;
    if ((adams = fopen("adams.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }


    for (uint i = 0; i < first_four.size(); ++i) {
        cout << first_four[i][0] << " " << first_four[i][1] << endl;
        s = "";
        s += to_string(first_four[i][0]) + "\n";
        fprintf(adams, s.c_str());
        s = "";
        s += to_string(first_four[i][1]) + "\n";
        fprintf(adams, s.c_str());
    }

    cout << "Погрешность: " << RungeRombergMethodAdams(first_four, h, 2) << endl;

    fclose(adams);
}
