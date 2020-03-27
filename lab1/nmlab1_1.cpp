
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;

vector<double> proiz(vector<vector<double>> matrixP, vector<double> matrixB, int m)
{
    vector<double> newB(m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            newB[i] += matrixP[i][j] * matrixB[j];
    return newB;
}

void decomLUP(vector<vector<double>> matrixA, vector<double> matrixB, vector<vector<double>> matrixU, vector<vector<double>> matrixL, vector<vector<double>> matrixP, int m, ofstream& output)
{
    int count = m;// var for "mu"
    for (int i = 0; i < m; i++)
    {
        //search of lead string
        double max = fabs(matrixU[i][i]);
        int row1 = i;
        int row2 = i;
        for (int j = m - 1; j >= i; j--)
            if (max < fabs(matrixU[j][i]))
            {
                max = matrixU[j][i];
                row2 = j;
            }
        //if we find new lead string we swap them in U & P
        if (row1 != row2)
        {
            for (int j = 0; j < m; j++)
            {
                double temp1 = matrixU[row1][j];
                matrixU[row1][j] = matrixU[row2][j];
                matrixU[row2][j] = temp1;

                double temp2 = matrixP[row1][j];
                matrixP[row1][j] = matrixP[row2][j];
                matrixP[row2][j] = temp2;
                if (j < i)
                {
                    double temp3 = matrixL[row1][j];
                    matrixL[row1][j] = matrixL[row2][j];
                    matrixL[row2][j] = temp3;
                }
            }
        }
        //main algorithm Gauss
        if (i != m - 1)
        {
            count--;
            vector<double> mu(count);
            int str = 1;
            if (i > 0) str += i;
            for (int j = 0; j < count; j++)
            {
                mu[j] = -(double)matrixU[j + str][i] / matrixU[i][i];
                matrixL[j + str][i] = -mu[j];
            }


            for (int j = 0; j < count; j++)
                for (int k = 0; k < m; k++)
                    matrixU[int(j + str)][k] += matrixU[i][k] * mu[j];
        }
    }
    //solving Lz=B
    matrixB = proiz(matrixP, matrixB, m);
    vector<double> z(m);
    z[0] = matrixB[0];
    for (int i = 1; i < m; i++)
    {
        double sum = 0;
        for (int j = 0; j < m; j++)
            sum += matrixL[i][j] * z[j];
        z[i] = matrixB[i] - sum;
    }
    //solving Ux = z
    vector<double> x(m);
    x[matrixA.size() - 1] = (double)z[matrixA.size() - 1] / matrixU[matrixA.size() - 1][matrixA.size() - 1];
    for (int i = m - 2; i >= 0; i--)
    {
        double sum = 0;
        for (int j = 0; j < m; j++)
            sum += matrixU[i][j] * x[j];
        x[i] = (double)1 / matrixU[i][i] * (z[i] - sum);
    }
    //output all data to the file output
    output << "Матрица перестановок P:\n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
            output << matrixP[i][j] << " ";
        output << "\n\n";
    }

    output << "Матрица L:\n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            output.width(15);
            output << matrixL[i][j] << " ";
        }
        output << "\n\n";
    }

    output << "Матрица U:\n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            output.width(15);
            output << matrixU[i][j] << " ";
        }
        output << "\n\n";
    }

    output << "Решение системы СЛАУ\n";
    for (int i = 0; i < m; i++)
    {
        output << "x" << i + 1 << " = ";
        output.width(2);
        output << x[i] << "\n";
    }
}

int main()
{
    setlocale(0, "Russian");
    int m = 4;
    vector<vector<double>> matrixA, matrixU, matrixL, matrixP;
    matrixA.assign(m, vector<double>(m));
    matrixU.assign(m, vector<double>(m));
    matrixL.assign(m, vector<double>(m));
    matrixP.assign(m, vector<double>(m));
    vector<double> matrixB(m);
    ifstream input("C:\\Users\\ksenia\\Desktop\\inp1_1.txt");
    ofstream output("C:\\Users\\ksenia\\Desktop\\outp1_1.txt");
    if (input.is_open())
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                input >> matrixA[i][j];

            }
            input >> matrixB[i];
        }
        matrixU = matrixA;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                matrixP[i][j] = 0;
                matrixL[i][j] = 0;
                if (i == j)
                {
                    matrixP[i][j] = 1;
                    matrixL[i][j] = 1;
                }
            }
        }
        output << "Матрица СЛАУ:\n";
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                //output << cout.width(3);
                output << matrixA[i][j] << " ";
            }
            output << "\n\n";
        }

        output << "Вектор значений:\n";
        for (int i = 0; i < m; i++)
            output << matrixB[i] << " ";

        output << "\n\nLUP-разложение матрицы А:\n\n";
        decomLUP(matrixA, matrixB, matrixU, matrixL, matrixP, m, output);

        output.close();
        input.close();
    }
    else
    {
        cout << "This file doesn't exist\n";
        cout << "Please check the name and way of file\n";
    }
}

