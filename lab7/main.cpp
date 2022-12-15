#include <iostream>
#include <vector>
#include <math.h>
#include <unistd.h>

using namespace std;

// x = [0, l]
// y = [0, k]
// hx - шаг по  x
// hy - шаг по y
// n - количество точек по x
// k - количество точек по y


const double PI = 3.141592653589793;

double f(double x, double y) {
    return 0;
}

void begin(vector<vector<double>> &u, int n, int k, double hx, double hy){
    for (int i = 0; i < n; ++i) {
        u[i][0] = exp(- i * hx) * cos(i * hx);
    }

    for (int i = 0; i < k; ++i) {
        u[0][i] = exp(- i * hy) * cos(i * hy);
    }    

    for (int i = 0; i < n; ++i) {
        u[i][k - 1] = 0;
    }

    for (int i = 0; i < k; ++i) {
        u[n - 1][i] = 0;   
    }
}

// void begin(vector<vector<double>> &u, int n, int k, double hx, double hy){
//     for (int i = 0; i < n; ++i) {
//         u[i][0] = exp(- i * hx) * cos(i * hx);
//     }

//     for (int i = 0; i < k; ++i) {
//         u[0][i] = cos(i * hy);
//     }    

//     for (int i = 0; i < n; ++i) {
//         u[i][k - 1] = 0;
//     }

//     for (int i = 0; i < k; ++i) {
//         u[n - 1][i] = 0;   
//     }

// }

vector<vector<double>> simpleIterations(vector<vector<double>> u, int n, int k, double hx, double hy, double a, double b1, double b2, double c, double eps, double relax, bool up_relax, int &p) {

    vector<vector<double>> nextu(n, vector<double>(k));

    double max = 1;

    while(max > eps) {

        max = -1;

        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < k - 1; ++j) {


                nextu[i][j] = (a * ((u[i - 1][j] + u[i + 1][j]) / (hx * hx) + (u[i][j - 1] + u[i][j + 1]) / (hy * hy)) -  
                              b1 * (u[i + 1][j] - u[i - 1][j]) / (2 * hx) - 
                              b2 * (u[i][j + 1] - u[i][j - 1]) / (2 * hy) - f(i * hx, j * hy)) / 
                              ( 2 * a * (1 / (hx * hx) + 1 / (hy * hy)) + c );

                // if (up_relax) {
                //     nextu[i][j] = u[i][j] + relax * (nextu[i][j] - u[i][j]); 
                // }

                double r = abs(u[i][j] - nextu[i][j]);

                if (max < r) {
                    max = r;
                }
            }   
        }  

        p += 1;
        swap(u, nextu);
        begin(u, n, k, hx, hy);
    }



    return u;
}

vector<vector<double>> seidel(vector<vector<double>> u, int n, int k, double hx, double hy, double a, double b1, double b2, double c, double eps, double relax, bool up_relax, int &p) {

    double max = 1;

    while(max > eps) {

        max = -1;

        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < k - 1; ++j) {

                double r = u[i][j];
                u[i][j] = (a * ((u[i - 1][j] + u[i + 1][j]) / (hx * hx) + (u[i][j - 1] + u[i][j + 1]) / (hy * hy)) -  
                              b1 * (u[i + 1][j] - u[i - 1][j]) / (2 * hx) - 
                              b2 * (u[i][j + 1] - u[i][j - 1]) / (2 * hy) - f(i * hx, j * hy)) / 
                              ( 2 * a * (1 / (hx * hx) + 1 / (hy * hy)) + c );

                if (up_relax) {
                    u[i][j] = r + relax * (u[i][j] - r); 
                }

                r = abs(u[i][j] - r);

                if (max < r) {
                    max = r;
                }
            }
        }

        p += 1;
    }

    return u;
}






int main() {
    double a, b1, b2, c, lx, ly, hx, hy, relax, eps;

    int n, k;
    cin >> lx >> n >> ly >> k >> a >> b1 >> b2 >> c >> relax >> eps;

    hx = lx / (n - 1);
    hy = ly / (k - 1);

    vector<vector<double>> u(n, vector<double>(k));
    begin(u, n, k, hx, hy);

    vector<vector<double>> sol(n, vector<double>(k));
    FILE * mainData;

    if ((mainData = fopen("data.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    fprintf(mainData,"%s", (to_string(a) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(b1) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(b2) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(c) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(hx) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(hy) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(k) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(n) + "\n").c_str());


    int p = 0;
    sol = simpleIterations(u, n, k, hx, hy, a, b1, b2, c, eps, relax, false, p);

    fprintf(mainData, "%s", (to_string(p) + "\n").c_str());

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }

    p = 0;

    sol = seidel(u, n, k, hx, hy, a, b1, b2, c, eps, relax, false, p);

    fprintf(mainData, "%s", (to_string(p) + "\n").c_str());

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }

    p = 0;

    sol = seidel(u, n, k, hx, hy, a, b1, b2, c, eps, relax, true, p);

    fprintf(mainData, "%s", (to_string(p) + "\n").c_str());

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }


    fclose(mainData);
}