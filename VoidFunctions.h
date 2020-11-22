//
// Created by fobos5892 on 11/21/20.
//

#ifndef UNTITLED2_VOIDFUNCTIONS_H
#define UNTITLED2_VOIDFUNCTIONS_H

#include "Source.h"

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

void GridGenerator();

double eps(double, double);

double deps_x(double, double);
double deps_y(double, double);

double mu(double, double);

double dmu_x(double, double);
double dmu_y(double, double);

double depsmu_x(double, double);
double depsmu_y(double, double);

double kappa(double, double, double, double);

double bf(double, double, double, double);
double dbf_x(double, double, double, double);
double dbf_y(double, double, double, double);

class VoidFunctions {

};


#endif //UNTITLED2_VOIDFUNCTIONS_H
