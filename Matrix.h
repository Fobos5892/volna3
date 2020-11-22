//
// Created by fobos5892 on 11/22/20.
//

#ifndef UNTITLED2_MATRIX_H
#define UNTITLED2_MATRIX_H

#include "VoidFunctions.h"

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

double matrix_ee(int, int, double, double, double, double, double, double**, int);
double matrix_mm(int, int, double, double, double, double, double, double**, int);
double matrix_em(int, int, double, double, double, double, double, double**, double**, int);
double matrix_me(int, int, double, double, double, double, double, double**, double**, int);

double matrix_solver();

void field_solver();

void Gauss(double**);

double Determinant(double**);

#endif //UNTITLED2_MATRIX_H
