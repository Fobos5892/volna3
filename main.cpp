#include "Source.h"
#include "VoidFunctions.h"
#include <linux/limits.h>
#include <string.h>
#include <time.h>
#include "Matrix.h"

double a = 1.0;
double b = a / 8;

double _gamma = 0;
double _omega = 0;

int m = 5;

double h = a / m;
double h_int = 10;

int p1 = 0;
int p2 = 0;

int N = p1 + p2;

vector<coordinates> GridE = vector<coordinates>();
vector<coordinates> GridM = vector<coordinates>();

char buffer_m[100];
char buffer_e[100];

int main() {
    std::cout << "Hello, welocme to a program!" << std::endl;

    char buffer[80];
    strcpy(buffer, "");
    char* buff2 = buffer;
    time_t rawtime;
    time ( &rawtime );
    struct tm * timeinfo;
    timeinfo = localtime ( &rawtime );
    strftime(buff2, 80, "%d%m%y", timeinfo);

    strcpy(buffer_m, "");
    strcat(buffer_m, buff2);
    strcat(buffer_m, "_fm");
    strcpy(buffer_e, "");
    strcat(buffer_e, buff2);
    strcat(buffer_e, "_fe");

    GridGenerator();

    N = p1 + p2;

    cout<<GridE[350].x<<endl;

    std::cout <<"Step of calculation " << h <<  endl;
    std::cout <<"The number of base elements "<<p1+p2<<" = "<< p1 << " + " << p2 << endl;
    std::cout <<"The dimension of matrix " << N <<"x"<<N<< endl;
    std::cout << "------------------------------------------------------------------------" << endl;

    ofstream file;
    file.open("znaki.xls"); //file of signum of system
    ofstream myfile;
    myfile.open("g_result.csv");

    double s = 0;
    double j = 2;
    _omega = j / 100;

    double gamma_h = 0.001;

    double gamma_min = 1 * _omega;
    double gamma_max = 3 * _omega;

    int gamma_count = 0;

    _gamma = gamma_max;

    double de_old = matrix_solver();
    //field_solver(p1, p2, N, h, a, b, omega, gamma, GridE, GridM, h_int);

    cout << "omega = " << _omega << " | gamma_min = " << gamma_min << " | gamma_max =  " << gamma_max << " | gamma_h = " << gamma_h << endl;

    for (double j = 2; j <= 100; j += 2)
    {

        for (double i = gamma_min; i <= gamma_max; i = i + gamma_h)
        {
            _gamma = gamma_max - i;

            double de_new = matrix_solver();
            field_solver();
            double count = de_new*de_old;
            file << de_new << "\t" << de_old << endl;
            if (count < 0)
            {
                field_solver();
                if (gamma_count > 0) break;
                gamma_count++;

                cout << _omega << " " << _gamma << " ----- " << s << endl;
                myfile << _omega << "," << _gamma << endl;

            }
            de_old = de_new;

        }
        cout << "------------------------------------------------------------------------" << endl;
    }

    std::cout << "------------------------------------------------------------------------" << endl;


    file.close();
    myfile.close();

    return 0;
}
