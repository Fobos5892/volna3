//
// Created by fobos5892 on 11/22/20.
//

#include "Matrix.h"
#include "Excel.h"
#include <cmath>

double matrix_ee(int ie, int je)
{
    double hx = h / h_int;
    double hy = h / h_int;

    double x10 = GridE[ie].x;
    double y10 = GridE[ie].y;

    double x20 = GridE[je].x;
    double y20 = GridE[je].y;

    double I1 = 0;

    if (sqrt((x10 - x20) * (x10 - x20) + (y10 - y20) * (y10 - y20)) <= 2 * h)
    {
        for (double x = x10 - h; x < x10 + h; x += hx)
        {
            for (double y = y10 - h; y < y10 + h; y += hy)
            {
                I1 += pow(_gamma, 2.0)
                        * eps(x + hx / 2.0, y + hy / 2.0)
                        * bf(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * bf(x + hx / 2.0, y + hy / 2.0, x20, y20)

                        + eps(x + hx / 2.0, y + hy / 2.0)
                        * dbf_x(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * dbf_x(x + hx / 2.0, y + hy / 2.0, x20, y20)
                        + dbf_y(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * dbf_y(x + hx / 2.0, y + hy / 2.0, x20, y20)

                        - pow(_omega, 2.0)
                        * pow(eps(x + hx / 2.0, y + hy / 2.0), 2.0)
                        * mu(x + hx / 2.0, y + hy / 2.0)
                        * bf(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * bf(x + hx / 2.0, y + hy / 2.0, x20, y20)

                        + (pow(_omega, 2.0)
                        * eps(x + hx / 2.0, y + hy / 2.0)
                        * bf(x + hx / 2.0, y + hy / 2.0, x10, y10))
                        / (pow(_gamma, 2.0)
                        - pow(_omega, 2.0)
                        * eps(x + hx / 2.0, y + hy / 2.0)
                        * mu(x + hx / 2.0, y + hy / 2.0))
                        * (depsmu_x(x + hx / 2.0, y + hy / 2.0)
                            * dbf_x(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        + depsmu_y(x + hx / 2.0, y + hy / 2.0)
                          * dbf_y(x + hx / 2.0, y + hy / 2.0, x10, y10));

            }
        }
    }

    return I1 * hx * hy;
}
double matrix_em(int ie, int jm)
{
    double hx = h / h_int;
    double hy = h / h_int;

    double x10 = GridE[ie].x;
    double y10 = GridE[ie].y;
    double x20 = GridM[jm].x;
    double y20 = GridM[jm].y;

    double I1 = 0;

    if (sqrt((x10 - x20) * (x10 - x20) + (y10 - y20) * (y10 - y20)) <= 2 * h)
    {
        for (double x = x10 - h; x < x10 + h; x += hx)
        {
            for (double y = y10 - h; y < y10 + h; y += hy)
            {
                I1 += _gamma * _omega
                        / (pow(_gamma, 2.0)
                        - pow(_omega, 2.0)
                        * eps(x + hx / 2.0, y + hy / 2.0)
                        * mu(x + hx / 2.0, y + hy / 2.0))
                        * bf(x + hx / 2.0, y + hy / 2.0, x20, y20)
                          * (depsmu_x(x + hx / 2, y + hy / 2)
                             * dbf_y(x + hx / 2, y + hy / 2, x10, y10)
                             - depsmu_y(x + hx / 2, y + hy / 2)
                               * dbf_x(x + hx / 2, y + hy / 2, x10, y10));
            }
        }
    }

    return I1 * hx * hy;
}
double matrix_me(int im, int je)
{
    double hx = h / h_int;
    double hy = h / h_int;

    double x10 = GridM[im].x;
    double y10 = GridM[im].y;
    double x20 = GridE[je].x;
    double y20 = GridE[je].y;

    double I1 = 0;

    if (sqrt((x10 - x20) * (x10 - x20) + (y10 - y20) * (y10 - y20)) <= 2 * h)
    {
        for (double x = x10 - h; x < x10 + h; x += hx)
        {
            for (double y = y10 - h; y < y10 + h; y += hy)
            {
                I1 += _gamma * _omega
                      / (pow(_gamma, 2.0)
                         - pow(_omega, 2.0)
                           * eps(x + hx / 2.0, y + hy / 2.0)
                           * mu(x + hx / 2.0, y + hy / 2.0))
                      * bf(x + hx / 2.0, y + hy / 2.0, x20, y20)
                      * (depsmu_x(x + hx / 2, y + hy / 2)
                         * dbf_y(x + hx / 2, y + hy / 2, x10, y10)
                         - depsmu_y(x + hx / 2, y + hy / 2)
                           * dbf_x(x + hx / 2, y + hy / 2, x10, y10));
            }
        }
    }

    return I1 * hx * hy;
}
double matrix_mm(int im, int jm)
{
    double hx = h / h_int;
    double hy = h / h_int;

    double x10 = GridM[im].x;
    double y10 = GridM[im].y;

    double x20 = GridM[jm].x;
    double y20 = GridM[jm].y;

    double I1 = 0;

    if (sqrt((x10 - x20) * (x10 - x20) + (y10 - y20) * (y10 - y20)) <= 2 * h)
    {
        for (double x = x10 - h; x < x10 + h; x += hx)
        {
            for (double y = y10 - h; y < y10 + h; y += hy)
            {
                I1 += pow(_gamma, 2.0)
                      * mu(x + hx / 2.0, y + hy / 2.0)
                      * bf(x + hx / 2.0, y + hy / 2.0, x10, y10)
                      * bf(x + hx / 2.0, y + hy / 2.0, x20, y20)

                      + mu(x + hx / 2.0, y + hy / 2.0)
                        * dbf_x(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * dbf_x(x + hx / 2.0, y + hy / 2.0, x20, y20)
                      + dbf_y(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * dbf_y(x + hx / 2.0, y + hy / 2.0, x20, y20)

                      - pow(_omega, 2.0)
                        * eps(x + hx / 2.0, y + hy / 2.0)
                        * pow(mu(x + hx / 2.0, y + hy / 2.0), 2.0)
                        * bf(x + hx / 2.0, y + hy / 2.0, x10, y10)
                        * bf(x + hx / 2.0, y + hy / 2.0, x20, y20)

                      + (pow(_omega, 2.0)
                         * eps(x + hx / 2.0, y + hy / 2.0)
                         * bf(x + hx / 2.0, y + hy / 2.0, x20, y20))
                        / (pow(_gamma, 2.0)
                           - pow(_omega, 2.0)
                             * eps(x + hx / 2.0, y + hy / 2.0)
                             * mu(x + hx / 2.0, y + hy / 2.0))
                        * (depsmu_x(x + hx / 2.0, y + hy / 2.0)
                           * dbf_x(x + hx / 2.0, y + hy / 2.0, x20, y20)
                           + depsmu_y(x + hx / 2.0, y + hy / 2.0)
                             * dbf_y(x + hx / 2.0, y + hy / 2.0, x10, y10));

            }
        }
    }

    return I1 * hx * hy;
}

double Determinant(double** matrix)
{
    int i, j, k;		/* Matrix subscripts */

    double tmp;	/* Temporary variable */

    for (i = 0; i < N; i++)
    {
        tmp = matrix[i][i];
        for (j = N; j >= i; j--)
            matrix[i][j] /= tmp;
        for (j = i + 1; j < N; j++)
        {
            tmp = matrix[j][i];
            for (k = N; k >= i; k--)
                matrix[j][k] -= tmp * matrix[i][k];
        }
    }

    for (i = 0; i < N; i++)
    {
        tmp *= matrix[i][i];
    }

    return tmp;
}

void Gauss(double** matrix)
{
    double* x, max;
    double* y = new double[N];

    for (int i = 0; i < N - 1; i++)
    {
        y[i] = 0.0;
    }
    y[N-1] = 1.0;

    int k, index;
    const double eps = 0.00001;  // точность
    x = new double[N];
    k = 0;
    while (k < N)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(matrix[k][k]);
        index = k;
        for (int i = k + 1; i < N; i++)
        {
            if (abs(matrix[i][k]) > max)
            {
                max = abs(matrix[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return;
        }
        for (int j = 0; j < N; j++)
        {
            double temp = matrix[k][j];
            matrix[k][j] = matrix[index][j];
            matrix[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < N; i++)
        {
            double temp = matrix[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < N; j++)
                matrix[i][j] = matrix[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < N; j++)
                matrix[i][j] = matrix[i][j] - matrix[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    saver_excel("Matrix_GAUSS_.xls", matrix);
    // обратная подстановка
    for (k = N - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - matrix[i][k] * x[k];
    }

    for (int i = 0; i < N; i++)
    {
        matrix[i][i] = x[i];
    }

    saver_excel("Matrix_GAUSS_inverse.xls", matrix);
}

double matrix_solver()
{
    double** Matrix = new double* [N];

    for (int i = 0; i < N; i++)
    {
        Matrix[i] = new double[N];
        for (int j = 0; j < N; j++)
        {
            Matrix[i][j] = 0;
        }
    }

    double res = 0;

    for (int i = 0; i < p1; i++)
    {
        for (int j = 0; j < p1; j++)
        {
            Matrix[i][j] = matrix_ee(i, j);
        }

        for (int j = 0; j < p2; j++)
        {
            Matrix[i][p1 + j] = matrix_em(i, j);

        }
    }

    for (int i = 0; i < p2; i++)
    {
        for (int j = 0; j < p1; j++)
        {
            Matrix[p1 + i][j] = matrix_me(i, j);
        }

        for (int j = 0; j < p2; j++)
        {
            Matrix[p1 + i][p1 + j] = matrix_mm(i, j);
        }
    }

    saver_excel("Matrix_solve.xls", Matrix);

    ofstream file;
    file.open("reserveMatrix.txt");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            file << Matrix[i][j] << "\t";
        }
        file << endl;
    }
    file.close();

    res = Determinant(Matrix);

    saver_excel("Matrix_Det.xls", Matrix);
    for (int i = 0; i < N; i++) { delete[]Matrix[i]; }

    delete[]Matrix;

    return res;
}

void field_solver()
{
    double** Matrix = new double* [N];

    for (int i = 0; i < N; i++)
    {
        Matrix[i] = new double[N];
    }

    ifstream file;
    file.open("reserveMatrix.txt");


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            file >> Matrix[i][j];
        }
    }

    file.close();

    saver_excel("Matrix_reading.xls", Matrix);

    Print_Basis(0, GridE);

    Print_Basis(1, GridM);

    Gauss(Matrix);

    Print_Field(0, GridE, Matrix);

    Print_Field(1, GridM, Matrix);



    for (int i = 0; i < N; i++) { delete[]Matrix[i]; }

    delete[]Matrix;
}