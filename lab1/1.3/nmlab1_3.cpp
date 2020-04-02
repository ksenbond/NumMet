
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int find(vector<int> row, int i, int m)
{
	int index = i;
	for (int j = 0; j < m; j++)
		if (i == row[j])
			index = j;
	return index;
}

void transpositionMatrix(vector<vector<double>> matrix, int m)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
		{
			matrix[j][i] = matrix[i][j];
			if (i > j)
				matrix[j][i] = 0;
		}
}

vector<vector<double>> get_minor(vector<vector<double>> newBE, vector<vector<double>> tempMatrix, int m, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < m; i++)
	{
		if (i != indRow)
		{
			for (int j = 0, kj = 0; j < m; j++)
			{
				if (j != indCol)
				{
					tempMatrix[ki][kj] = newBE[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
	return tempMatrix;
}

double determinate(vector<vector<double>> newBE, int m)
{
	double det = 0;
	int k = 1;    //степень  
	if (m == 1)
		det = newBE[0][0];
	if (m == 2)
		det = newBE[0][0] * newBE[1][1] - newBE[1][0] * newBE[0][1];
	else
	{
		for (int i = 0; i < m; i++)
		{
			int n = m - 1;
			vector<vector<double>> tempMatrix;
			tempMatrix.assign(n, vector<double>(n));
			for (int j = 0; j < n; j++)
				for (int d = 0; d < n; d++)
					tempMatrix[j][d] = 0;
			tempMatrix = get_minor(newBE, tempMatrix, m, 0, i);
			det += k * newBE[0][i] * determinate(tempMatrix, n);
			k = -k;
			tempMatrix.clear();
		}
	}
	return det;
}

vector<vector<double>> inverse(vector<vector<double>> newBE, double det, int m)
{
	vector<vector<double>> inverse;
	inverse.assign(m, vector<double>(m));
	for (int i = 0; i < m; i++) 
	{
		for (int j = 0; j < m; j++) 
		{
			int n = m - 1;
			vector<vector<double>> tempMatrix;
			tempMatrix.assign(n, vector<double>(n));
			tempMatrix = get_minor(newBE, tempMatrix, m, i, j);
			inverse[i][j] = pow(-1.0, i + j + 2) * determinate(tempMatrix, n) / det;
			tempMatrix.clear();
		}
	}
	transpositionMatrix(inverse, m);
	return inverse;
}

void zeydel(vector<vector<double>> matrixA, vector<double> matrixB, int m, double eps, fstream& out)
{
	vector<vector<double>> alpha;
	alpha.assign(m, vector<double>(m));
	vector<double> beta(m);
	vector<vector<double>> matrixC;
	matrixC.assign(m, vector<double>(m));
	vector<vector<double>> newBE;
	newBE.assign(m, vector<double>(m));
	//creation of matrix alpha &  * (E-Bbeta from book
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++) 
		{
			if (i == j)
				alpha[i][j] = 0;
			else
				alpha[i][j] = -matrixA[i][j] / matrixA[i][i];
		}
		beta[i] = matrixB[i] / matrixA[i][i];
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i > j)
				newBE[i][j] = -alpha[i][j];
			if (i < j)
				matrixC[i][j] = alpha[i][j];
			if (i == j)
			{
				newBE[i][j] = 1;
				matrixC[i][j] = 0;
			}
		}
	}
	double det = determinate(newBE, m);
	vector<vector<double>> inverseBE;
	inverseBE.assign(m, vector<double>(m));
	if (det != 0)
		inverseBE = inverse(newBE, det, m); 
	vector<vector<double>> newBEC;
	newBEC.assign(m, vector<double>(m));
	vector<double> newBeta(m);
	for (int i = 0; i < m; i++) 
	{
		for (int j = 0; j < m; j++) 
		{
			newBeta[i] += inverseBE[i][j] * beta[j];
			newBEC[i][j] = 0;
			for (int k = 0; k < m; k++) 
			{
				newBEC[i][j] += inverseBE[i][k] * matrixC[k][j];
			}
		}
	}
	int flag = 0;
	int counter = 0;
	vector<double> x(m);
	while (flag == 0) 
	{
		for (int i = 0; i < m; i++)
		{
			if (i == 0 && counter == 0)
				for (int j = 0; j < m; j++)
					x[j] = newBeta[j];
			else
			{
				double sum = 0;
				for (int j = 0; j < m; j++)
					sum += newBEC[i][j] * x[j];
				double temp = x[i];
				x[i] = newBeta[i] + sum;
				if (abs(x[i] - temp) < eps)
					flag = 1;
			}
		}
		counter++;
	}
	for (int i = 0; i < m; i++)
	{
		out << "x" << i + 1 << " = ";
		out.width(5);
		out << x[i] << "\n";
	}
	out << "amount of iterations = " << counter << "\n";
}

int main()
{
	fstream in("input.txt");
	fstream out("output.txt");
	double eps;
	int m;
	if (in.is_open())
	{
		in >> eps >> m;
		vector<vector<double>> matrixA;
		matrixA.assign(m, vector<double>(m));
		vector<double> matrixB(m);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
				in >> matrixA[i][j];
			in >> matrixB[i];
		}
		//max a[i][i]
		vector<int> row(m);
		for (int i = 0; i < m; i++)
		{
			double maxi = 0;
			for (int j = 0; j < m; j++)
			{
				if (abs(matrixA[i][j]) > maxi)
				{
					maxi = abs(matrixA[i][j]);
					row[i] = j;
				}
			}
		}
		//diagonal prevalence
		for (int i = 0; i < m; i++)
		{
			if (i != row[i])
			{
				int index = find(row, i, m);
				for (int j = 0; j < m; j++)
				{
					double temp1 = matrixA[i][j];
					matrixA[i][j] = matrixA[index][j];
					matrixA[index][j] = temp1;
				}
				double temp2 = matrixB[i];
				matrixB[i] = matrixB[index];
				matrixB[index] = temp2;

				int temp3 = row[i];
				row[i] = row[index];
				row[index] = temp3;
			}
		}
		//main algorithm
		zeydel(matrixA, matrixB, m, eps, out);
		in.close();
		out.close();
	}
	else
	{
		cout << "Such file doesn't exist";
	}
}