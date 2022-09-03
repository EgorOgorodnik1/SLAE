#include <iostream>
#include <cmath>

using namespace std;

const double eps = 0.001;

const unsigned int len = 3;

void matrix_output(double matrix[][len], double* b)
{
	for (int i = 0; i < len; i++){
		for (int j = 0; j < len; j++)
			cout << matrix[i][j] << "\t";
		cout << "|" << b[i] << endl;
	}
	cout << endl;
}

void matrixTranspose(double matrix[][len], double matrix_TP[][len])
{
	for (int i = 0; i < len; i++){
		for (int j = 0; j < len; j++)
			matrix_TP[i][j] = matrix[i][j];
	}

	double temp;
	for (int i = 0; i < len; ++i){
		for (int j = i; j < len; ++j){
			temp = matrix_TP[i][j];
			matrix_TP[i][j] = matrix_TP[j][i];
			matrix_TP[j][i] = temp;
		}
	}
}

void matrixTransformation(double matrix[][len], double matrix_TP[][len], double matrix_TF[][len], double* b_TF, double* b)
{
	for (int i = 0; i < len; i++){
		
		b_TF[i] = 0;
		for (int j = 0; j < len; j++){
			matrix_TF[i][j] = 0;

			b_TF[i] += matrix_TP[i][j] * b[j];

			for (int k = 0; k < len; k++){
				matrix_TF[i][j] += matrix_TP[i][k] * matrix[k][j];
			}
		}
	}
	matrix_output(matrix_TF, b_TF);
}

void methodSeidel(double matrix[][len], double* b, double matrix_TF[][len], double* b_TF, double* x_next, double* x_pr, double p)
{
	int counter = 0;

	double delta_X = 0;

	double discrepancy_norm = 0;

	double first_sum, second_sum, delta_sum, discrepancy_sum;

	do
	{

		for (int i = 0; i < len; i++){
			first_sum = 0;
			for (int j = 0; j < i; j++){
				first_sum += matrix_TF[i][j] * x_next[j];
			}

			second_sum = 0;
			for (int j = i + 1; j < len; j++){
				second_sum += matrix_TF[i][j] * x_pr[j];
			}
			x_next[i] = ((1 + p) * (b_TF[i] - first_sum - second_sum) / matrix_TF[i][i]) - p * x_pr[i];
		}

		delta_sum = 0;
		for (int k = 0; k < len; k++){
			delta_sum += (x_next[k] - x_pr[k]) * (x_next[k] - x_pr[k]);
		}
		delta_X = sqrt(delta_sum);

		discrepancy_sum = 0;
		for (int i = 0; i < len; i++){
			double discrepancy = 0;
			for (int j = 0; j < len; j++){
				discrepancy += matrix[i][j] * x_next[j];
			}
			discrepancy_sum += (discrepancy - b[i]) * (discrepancy - b[i]);
		}
		discrepancy_norm = sqrt(discrepancy_sum);

		for (int l = 0; l < len; l++)
			x_pr[l] = x_next[l];

		counter++;
	} 
	while (delta_X >= eps);

	cout << "||AX-B||" << "\t" << "||X[k+1]-X[k]|| < e" << endl << endl;
	cout << discrepancy_norm << "         " << delta_X << endl << endl;
	cout << "Число итераций - " << counter << endl << endl;
	cout << "Решение: X = " << '(';
	for (int i = 0; i < len; i++)
		cout << x_next[i] << ' ';
	cout << ')' << endl;
	for (int i = 0; i < 40; i++)
		cout << '-';
	cout << endl;
}

int main()
{
	setlocale(LC_ALL, "Rus");

	float p = 0;

	double matrix[][len] = {
		{8, 10, 0},
		{5, 1, -3},
		{0, 3,  9}
	};

	double b[len] = { 36,3,33};

	double matrix_TP[len][len];

	double matrix_TF[len][len];

	double b_TF[len];

	double x_next[len];

	double x_pr[len];

	cout << "Исходная СЛАУ"<<endl;

	matrix_output(matrix, b);

	matrixTranspose(matrix, matrix_TP);

	cout << "Приведение СЛАУ к симметричному виду" << endl;

	matrixTransformation(matrix, matrix_TP, matrix_TF, b_TF, b);

	cout << "Начальное приближение: X[0] = ( ";
	for (int i = 0; i < len; i++){
		x_pr[i] = 0;
		cout << x_pr[i] << ' ';
	}
	cout << ')' << endl;

	for (int i = 0; i < 40; i++)
		cout << '-';
	cout << endl << endl;

	while (p < 1)
	{
		cout << "Ускоряющий параметр p = " << p << ':' << endl << endl;

		for (int i = 0; i < len; i++)
			x_pr[i] = b[i] / matrix[i][i];

		for (int i = 0; i < len; i++)
			x_next[i] = 0;

		methodSeidel(matrix, b, matrix_TF, b_TF, x_next, x_pr, p);

		cout << endl;

		p += 0.1;
	}

	return 0;
}
