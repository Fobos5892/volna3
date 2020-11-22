//
// Created by fobos5892 on 11/22/20.
//

#ifndef UNTITLED2_EXCEL_H
#define UNTITLED2_EXCEL_H

#include "Source.h"

#define MAX_PATH 260

extern double a;
extern double b;

extern double _gamma;
extern double _omega;

extern int m;

extern double h;
extern double h_int;

extern int p1;
extern int p2;

extern int N;

extern vector<coordinates> GridE;
extern vector<coordinates> GridM;

void saver_excel(char*, double**);

void WriteVector(char*, char* , vector<coordinates>, double*);

void Print_Field(int, vector<coordinates>, double **);

void Print_Basis(int, vector<coordinates>);

#endif //UNTITLED2_EXCEL_H
