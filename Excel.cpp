//
// Created by fobos5892 on 11/22/20.
//

#include "Excel.h"
#include "VoidFunctions.h"
#include <cstdlib>
#include <iomanip>
#include <cstring>

extern char buffer_m[100];
extern char buffer_e[100];

void saver_excel(char* name, double** X)
{
    ofstream file;

    file.open(name);

    int l = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            file << X[i][j] << "\t";
        }
        file << "\n";
    }

    file.close();
}

void WriteVector(char* filname1, char* filename2, vector<coordinates> Grid, double* U)
{

    ofstream myfile1, myfile2;

    cout << sizeof(U) << endl;

    myfile1.open(filname1);
    myfile2.open(filename2);

    char* HHH = "\"\"\,\"\"\,\"\"\n";

    myfile1 << HHH;
    myfile2 << HHH;

    ofstream fileE, fileM;

    for (int i = 0; i < p1; i++)
    {
        myfile1 << setprecision(5) << Grid[i].x << "," << setprecision(5) << Grid[i].y << "," << setprecision(5) << U[i] << endl;
    }

    for (int i = p1; i < p1 + p2; i++)
    {
        myfile2 << setprecision(5) << Grid[i].x << "," << setprecision(5) << Grid[i].x << "," << setprecision(5) << U[i] << endl;
    }

    myfile1.close();
    myfile2.close();

}

void Print_Field(int param, vector<coordinates> Grid, double **Matrix)
{
    int start = 0;
    int finish = p1;
    char *filename;
    char *folder;

    strcpy(filename, "");
    if (param == 0)
    {
        filename = "\\fe_field_";
        folder = "_fe";
    }
    else
    {
        filename = "\\fe_field_";
        folder = "_fe";
    }

    char filename_const[260];
    char chh[40];
    ofstream myfile;

    strcpy(filename_const, "");
    strcat(filename_const, buffer_m);
    strcat(filename_const, filename);
    sprintf(chh, "%f", _gamma);
    strcat(filename_const, chh);
    strcat(filename_const, "_");
    sprintf(chh, "%f", _omega);
    strcat(filename_const, chh);
    strcat(filename_const, ".csv");

    myfile.open(filename_const);

    char* HHH = "\"\"\,\"\"\,\"\"\n";

    myfile << HHH;

    double hx = h / h_int; double hy = h / h_int;

    for (int i = start; i < finish; i++)
    {
        myfile << Grid[i].x << "," << Grid[i].y << "," << Matrix[i][i] << endl;
    }

    myfile.close();
}

void Print_Basis(int param, vector<coordinates> Grid)
{
    int start = 0;
    int finish = p1;
    char* filename;
    char* folder;
    switch (param)
    {
        case 0:
            filename = "\\fe_basis_";
            folder = "_fe";
            break;
        case 1:
            filename = "\\fm_basis_";
            folder = "_fm";
            start = 0;
            finish = p2;
        default:
            break;
    }

    char filename_const[260];
    char chh[40];
    ofstream myfile;

    strcpy(filename_const, "");
    strcat(filename_const, buffer_m);
    strcat(filename_const, filename);
    sprintf(chh, "%f", _gamma);
    strcat(filename_const, chh);
    strcat(filename_const, "_");
    sprintf(chh, "%f", _omega);
    strcat(filename_const, chh);
    strcat(filename_const, ".csv");

    myfile.open(filename_const);

    char* HHH = "\"\"\,\"\"\,\"\"\n";

    myfile << HHH;

    double hx = h / h_int;
    double hy = h / h_int;

    for (int i = start; i < finish; i++)
    {
        myfile << setprecision(5) << Grid[i].x << "," << setprecision(5)
            << Grid[i].y<< "," << setprecision(5)
            << bf(Grid[i].x, Grid[i].y, Grid[i].x, Grid[i].y) << endl;
    }

    myfile.close();
}
