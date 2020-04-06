
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

void findMax(vector<vector<double>> matrix, int size, int *row, int *col, double *max)
{
    *max = 0;
    vector<double> diag(size);
    for (int i = 0; i < size; i++)
        diag[i] = abs(matrix[i][i]);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (j > i && (abs(matrix[i][j]) > *max) && (abs(matrix[i][j]) != diag[i]) && (abs(matrix[i][j]) != diag[j]))
            {
                *max = abs(matrix[i][j]);
                *row = i;
                *col = j;
            }
}

double sumSquares(vector<vector<double>> matrix, int size)
{
    double sum = 0;
    for (int i = 0; i < size; i++)
        for (int j = i; j < size; j++)
            if (j > i)
                sum += matrix[i][j] * matrix[i][j];
    return sqrt(sum);

}
 
vector<vector<double>> matrixMultiplication(vector<vector<double>> matrixA, vector<vector<double>> matrixB, int size)
{
    vector<vector<double>> result;
    result.assign(size, vector<double>(size));
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
    return result;
}

void algorithmJakobi(vector<vector<double>> matrix, int size, double eps, fstream &out)
{
    int row = 0,
        col = 0,
        counter = 0,
        flag = 0;
    double max = 0,
           fi = 0;
    vector<vector<double>> matrixU;
    vector<vector<double>> transposU;
    vector<vector<double>> eigenvector;
    //vector<vector<double>> lambda;
    matrixU.assign(size, vector<double>(size));
    transposU.assign(size, vector<double>(size));
    eigenvector.assign(size, vector<double>(size));
    //lambda.assign(size, vector<double>(size));
    double sum = 10;
    do//main algorithm
    {
        findMax(matrix, size, &row, &col, &max); // max non-diagonal elem of matrix
        if (matrix[row][row] != matrix[col][col])
            fi = (double)1 / 2 * atan(2 * matrix[row][col] / (matrix[row][row] - matrix[col][col]));
        else
            fi = M_PI / 4;
        if (flag == 0)
        {
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    if (i == j)
                        eigenvector[i][j] = 1;
                    else
                        eigenvector[i][j] = 0;
                }
        }

        flag = 1;
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
            {
                matrixU[i][j] = 0;
                transposU[j][i] = 0;
                if (i != row && j != col && i == j)
                {
                    matrixU[i][j] = 1;
                    transposU[j][i] = 1;
                }
                if (i == row && j == col)
                {
                    matrixU[i][j] = -sin(fi);
                    transposU[j][i] = -sin(fi);
                }
                if (i == col && j == row)
                {
                    matrixU[i][j] = sin(fi);
                    transposU[j][i] = sin(fi);
                }
                if ((i == row && j == row) || (i == col && j == col))
                {
                    matrixU[i][j] = cos(fi);
                    transposU[j][i] = cos(fi);
                }
            }
        eigenvector = matrixMultiplication(eigenvector, matrixU, size);
        matrix = matrixMultiplication(matrixMultiplication(transposU, matrix, size), matrixU, size); // L = 
        counter++;
        sum = sumSquares(matrix, size);
    } while (sum > eps);
    vector<double> lambda(size);
    for (int i = 0; i < size; i++)
    {
        lambda[i] = matrix[i][i];
        out << "l" << i+1 << " = " << lambda[i] << "\n";
    }
    for (int i = 0; i < size; i++)
    {
        out << i + 1 << "й собственный вектор имеет вид:\n";
        for (int j = 0; j < size; j++)
        {
            out.width(15);
            out << eigenvector[j][i] << "\n";
        }
        out << "\n\n";
    }
    out << "количество итераций = " << counter;
}

int main()
{
    setlocale(0, "Russian");
    vector<vector<double>> matrixA;
    fstream in("in.txt");
    fstream out("out.txt");
    double eps;
    int size;
    if (in.is_open())
    {
        in >> eps >> size;
        matrixA.assign(size, vector<double>(size));
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                in >> matrixA[i][j];
        algorithmJakobi(matrixA, size, eps, out);
        in.close();
        out.close();
    }
    else
    {
        cout << "No such file";
    }
}
