
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

int sign(double var)
{
    if (var > 0) return 1;
    if (var < 0) return -1;
    if (var == 0) return 0;
}

vector<vector<double>> getMatrixH(vector<double> v, int size)
{
    vector<vector<double>> result;
    result.assign(size, vector<double>(size));
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += v[i] * v[i];
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i != j)
                result[i][j] = -(double)2/sum * v[i] * v[j];
            else 
                result[i][j] = 1 - (double)2/sum * v[i] * v[j];
        }
    }
    return result;
}

vector<vector<double>> matrixMultiplication(vector<vector<double>> matrixA, vector<vector<double>> matrixB, int size)
{
    vector<vector<double>> result;
    result.assign(size, vector<double>(size));
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            for (int k = 0; k < size; k++)
                result[i][j] += matrixA[i][k] * matrixB[k][j];
    return result;
}

vector<vector<double>> decompostionQR(vector<vector<double>> matrix, int size)
{
    vector<vector<double>> matrixR;
    vector<vector<double>> matrixH;
    vector<vector<double>> matrixQ;
    matrixR.assign(size, vector<double>(size));
    matrixH.assign(size, vector<double>(size));
    matrixQ.assign(size, vector<double>(size));
    int counter = 0;
    matrixR = matrix;
    for (int i = 0; i < size; i++)
    {
        vector<double> v(size);
        for (int j = 0; j < size; j++)
        {
            if (j > counter)
            {
                v[j] = matrixR[j][i];
            }
            if (j < counter)
            {
                v[j] = 0;
            }
            if (j == counter)
            {
                double sum = 0;
                for (int k = counter; k < size; k++)
                    sum += matrixR[k][i] * matrixR[k][i];
                v[j] = matrixR[i][i] + sign(matrixR[i][i]) * sqrt(sum);
            }
        }
        vector<vector<double>> matrixH;
        matrixH.assign(size, vector<double>(size));
        matrixH = getMatrixH(v, size);
        if (counter == 0)
            matrixQ = matrixH;
        else
            matrixQ = matrixMultiplication(matrixQ, matrixH, size);
        counter++;
        vector<vector<double>> matrixA;
        matrixA.assign(size, vector<double>(size));
        matrixA = matrixMultiplication(matrixH, matrixR, size);
        matrixR = matrixA;
    }
    vector<vector<double>> newA;
    newA.assign(size, vector<double>(size));
    newA = matrixMultiplication(matrixR, matrixQ, size);
    return newA;
}

void findlambda(vector<vector<double>> matrix, int size, double sum1, double sum2, double eps, fstream &out)
{
    //если все поддиагональные элементы равны 0
    if (sum1 == 0 && sum2 == 0)
    {
        vector<double> lambda(size);
        out << "Eigen values of matrix:\n";
        for (int i = 0; i < size; i++)
        {
            lambda[i] = matrix[i][i];
            out << "l" << i + 1 << " = " << lambda[i] << "\n";
        }
    }
    else
    {
        if (sqrt(sum1) < eps)
        {
            out << "Eigen values of matrix:\n";
            double lambda = matrix[0][0];
            out << "l1 = " <<  lambda << "\n";
            int a;
            if (sign(matrix[1][1]) != sign(matrix[2][2]))
                a = -1;
            else
                a = 1;
            double b = (matrix[1][1] + matrix[2][2]) * (matrix[1][1] + matrix[2][2]);
            double D = b - 4 * a * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
            out << "l2 = " << (double)-b / 2 / a << "+ i * " << sqrt(D) / 2 / a << "\n";
            out << "l3 = " << (double)-b / 2 / a << "- i * " << sqrt(D) / 2 / a << "\n";
        }
        else
        {
            out << "Eigen values of matrix:\n";
            double lambda = matrix[size - 1][size - 1];
            out << "l1 = " << lambda << "\n";
            int a;
            if (sign(matrix[1][1]) != sign(matrix[0][0]))
                a = -1;
            else
                a = 1;
            double b = (matrix[1][1] + matrix[0][0]) * (matrix[1][1] + matrix[0][0]);
            double D = b - 4 * a * (matrix[1][1] * matrix[0][0] - matrix[1][0] * matrix[0][1]);
            out << "l2 = " << (double)-sqrt(b) / 2 / a << " + i * " << sqrt(-D) / 2 / a << "\n";
            out << "l3 = " << (double)-sqrt(b) / 2 / a << " - i * " << sqrt(-D) / 2 / a << "\n";
        }
    }
    out.close();
}

void algorithm(vector<vector<double>> matrix, int size, double eps, fstream& out)
{
    int counter = 0,
        row = 0;
    vector<vector<double>> matrixA;
    matrixA.assign(size, vector<double>(size));
    matrixA = matrix;
    double sum1,
           sum2,
           sum;
    do
    {
        sum1 = 0;
        sum2 = 0;
        sum = 0;
        matrixA = decompostionQR(matrixA, size);
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
            {
                if (j > i && i == 0)
                    sum1 += matrixA[j][i] * matrixA[j][i];
                if (j > i && i == 1)
                    sum2 += matrixA[j][i] * matrixA[j][i];
            }
        sum = min(sqrt(sum1), sqrt(sum2));
        counter++;
    } while (sum > eps);
    vector<double> lambda(size);
    out << "Amount of iterations with eps = " << eps << " is " << counter << "\n";
    findlambda(matrixA, size, sum1, sum2, eps, out);
    out.close();
}

int main()
{
    fstream in("in.txt");
    fstream out("out.txt");
    double eps;
    int size;
    vector<vector<double>> matrix;
    if (in.is_open())
    {
        in >> eps >> size;
        matrix.assign(size, vector<double>(size));
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                in >> matrix[i][j];
        algorithm(matrix, size, eps, out);
        in.close();
        out.close();
    }
    else
    {
        out << "Such file doesn't exist";
    }
}