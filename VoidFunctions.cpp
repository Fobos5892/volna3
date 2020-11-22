//
// Created by fobos5892 on 11/21/20.
//

#include "VoidFunctions.h"


void GridGenerator()
{
    ofstream fileE, fileM;
    fileE.open("GridE.xls");
    fileM.open("GridM.xls");
    for (int i = 0; i <= 2 * m; i++)
    {
        for (int j = 0; j <= 2 * m; j++)
        {
            double x = -a + h * j;
            double y = -a + h * i;

            if ((abs(x) <= a) && (abs(y) <= a))
            {
                coordinates xy(x, y);
                GridM.push_back(xy);

                fileM << GridM[p1].x << "\t" << GridM[p1].y << endl;

                p1 = p1 + 1;
            }


            if ((abs(x) < a) && (abs(y) < a))
            {
                coordinates xy(x, y);
                GridE.push_back(xy);

                fileE << GridE[p1].x << "\t" << GridE[p2].y << endl;
                p2 = p2 + 1;
            }
        }


    }

    fileE.close();
    fileM.close();
}

double eps(double x, double y)
{
    return 4 + x * x * y * y;
}
double deps_x(double x, double y)
{
    return 2 * x * y * y;
}
double deps_y(double x, double y)
{
    return 2 * x * x * y;
}

double mu(double x, double y)
{
    return 4 + x * x * y * y;
}
double dmu_x(double x, double y)
{
    return 2 * x * y * y;
}
double dmu_y(double x, double y)
{
    return 2 * x * x * y;
}

double depsmu_x(double x, double y)
{
    return dmu_x(x, y) * eps(x, y) - deps_x(x, y) * mu(x, y);
}

double depsmu_y(double x, double y)
{
    return dmu_y(x, y) * eps(x, y) - deps_y(x, y) * mu(x, y);
}

double kappa(double x, double y)
{
    double res = _omega * _omega * eps(x, y) - _gamma * _gamma;
    return res;
}

double bf(double x0, double y0, double x, double y)
{
    double ax = x0 - h;
    double bx = x0;
    double cx = x0 + h;

    double ay = y0 - h;
    double by = y0;
    double cy = y0 + h;

    double x1 = (x - ax) / h;
    double x2 = (cx - x) / h;

    double y1 = (y - ay) / h;
    double y2 = (cy - y) / h;

    double fx, fy = 0;

    //basis x
    if (x < ax || x > cx || abs(x) > a)
    {
        fx = 0.0;
    }

    if ((x >= ax) && (x < bx))
    {
        fx = x1;
    }

    if (x == bx)
    {
        fx = 1.0;
    }

    if ((x > bx) && (x <= cx))
    {
        fx = x2;
    }

    //basis y
    if (y < ay || y > cy || abs(y) > a)
    {
        fy = 0.0;
    }

    if ((y >= ay) && (y < by))
    {
        fy = y1;
    }

    if (y == by)
    {
        fy = 1.0;
    }

    if ((y > by) && (y <= cy))
    {
        fy = y2;
    }

    return fx * fy;
}
double dbf_x(double x0, double y0, double x, double y)
{
    double ax = x0 - h;
    double bx = x0;
    double cx = x0 + h;

    double ay = y0 - h;
    double by = y0;
    double cy = y0 + h;

    double x1 = 1.0 / h;
    double x2 = - 1.0 / h;

    double y1 = 1.0 / h;
    double y2 = - 1.0 / h;

    double fx, fy = 0;

    //basis x
    if (x < ax || x > cx || abs(x) > a)
    {
        fx = 0.0;
    }

    if ((x >= ax) && (x < bx))
    {
        fx = x1;
    }

    if (x == bx)
    {
        fx = 1.0;
    }

    if ((x > bx) && (x <= cx))
    {
        fx = x2;
    }

    //basis y
    if (y < ay || y > cy || abs(y) > a)
    {
        fy = 0.0;
    }

    if ((y >= ay) && (y < by))
    {
        fy = y1;
    }

    if (y == by)
    {
        fy = 1.0;
    }

    if ((y > by) && (y <= cy))
    {
        fy = y2;
    }

    return fx * fy;
}
double dbf_y(double x0, double y0, double x, double y)
{
    double res = 0;

    double ax = x0 - h;
    double bx = x0;
    double cx = x0 + h;

    double ay = y0 - h;
    double by = y0;
    double cy = y0 + h;

    double x1 = x / (bx - ax) + ax / (ax - bx);
    double x2 = x / (bx - cx) + cx / (cx - bx);

    double y1 = 1 / (by - ay);
    double y2 = 1 / (by - cy);

    double fx, fy = 0;

    if (x < ax)
    {
        fx = 0.0;
    }

    if ((x >= ax) && (x < bx))
    {
        fx = x1;
    }

    if (x == bx)
    {
        fx = 1.0;
    }

    if ((x > bx) && (x <= cx))
    {
        fx = x2;
    }

    if (x > cx)
    {
        fx = 0.0;
    }


    if (y < ay)
    {
        fy = 0.0;
    }

    if ((y >= ay) && (y < by))
    {
        fy = y1;
    }

    if (y == by)
    {
        fy = 1.0;
    }

    if ((y > by) && (y <= cy))
    {
        fy = y2;
    }

    if (y > cy)
    {
        fy = 0.0;
    }

    res = fx * fy;

    if ((abs(x) > a) || (abs(y) > a)) { res = 0; }
    //if ((abs(x) < b) && (abs(y) < b)) { res = 0; }

    return res;
}