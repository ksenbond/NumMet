#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//main run-through algorithm
void rtm(vector<vector<double>> matrixA, vector<double> matrixB, int m, fstream& out)
{
	vector<double> p(m);
	vector<double> q(m);
	vector<double> x(m);
	for (int i = 0; i < m; i++)
	{
		if (i == 0)
		{
			p[i] = (double) -matrixA[i][i + 1] / matrixA[i][i];
			q[i] = (double)matrixB[i] / matrixA[i][i];
		}

		if (i != m - 1 && i != 0) 
		{
			p[i] = (double) -matrixA[i][i + 1] / (matrixA[i][i] + matrixA[i][i - 1] * p[i - 1]);
			q[i] = (double)(matrixB[i] - matrixA[i][i - 1] * q[i - 1]) / (matrixA[i][i] + matrixA[i][i - 1] * p[i - 1]);
			
		}

		if (i == m - 1)
		{
			p[i] = 0;
			q[i] = (double)(matrixB[i] - matrixA[i][i - 1] * q[i - 1]) / (matrixA[i][i] + matrixA[i][i - 1] * p[i - 1]);
		}
	}
	x[m - 1] = q[m - 1];
	for (int i = m - 2; i >= 0; i--)
		x[i] = p[i] * x[i + 1] + q[i];
	for (int i = 0; i < m; i++)
	{
		out << "x" << i + 1 << "=";
		out.width(2);
		out << x[i] << "\n";
	}
}

int main()
{
	int m = 5;
	fstream in("input.txt");
	fstream out("output.txt");
	vector <vector <double>> matrixA;
	matrixA.assign(m, vector<double>(m));
	vector<double> matrixB(m);
	if (in.is_open())
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				in >> matrixA[i][j];
			}
			in >> matrixB[i];
		}
		out << "Матрица СЛАУ:\n";
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				out.width(2);
				out << matrixA[i][j] << " ";
			}
			out << "\n\n";
		}
		out << "\nВектор правых частей\n";
		for (int i = 0; i < m; i++)
			out << matrixB[i] << " ";
		out << "\n";
		rtm(matrixA, matrixB, m, out);
		in.close();
		out.close();
	}
	else
	{
		cout << "File 'input.txt' doesn't exist";
	}
}
