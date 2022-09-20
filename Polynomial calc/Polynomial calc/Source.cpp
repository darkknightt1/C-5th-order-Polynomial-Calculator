

#include <iostream>
#include <cmath>
#include <string>
#include<complex>








using namespace std;
complex <double>cbrt(complex <double> a, int n)
{
	const double TWOPI = 2.0 * 3.141592653589793238462643383279502884;

	double rho = pow(abs(a), 1.0 / 3.0);
	double theta = ((TWOPI * n) + arg(a)) / 3.0;
	return complex<double>(rho * cos(theta), rho * sin(theta));
}

double calculate(double arr[], int degree, double x)
{
	double sum = 0;
	double power = degree;
	for (int i = 0; i < degree; i++)
	{
		sum += arr[i] * pow(x, power);
		power--;
		if (power == 0)
		{
			sum += arr[i + 1];
			break;
		}
	}
	return sum;
}

double sum_int(char arr[], int start, int end)
{
	double sum = 0;
	int power = 0;

	for (int i = end; i >= start; i--)
	{
		if (arr[i] == ' ')
		{
			continue;
		}
		sum += (arr[i] - '0') * pow(10, power);
		power++;
	}
	/*if (sum == 0)

		sum = 1;

	else if (sum == -0)

		sum = -1;
*/
	return sum;
}

double decimal_sum(char arr[], int start, int end)
{
	double sum = 0;
	int power = -1;
	for (int i = start; i <= end; i++)
	{
		if (arr[i] == ' ')
		{
			continue;
		}
		sum += (arr[i] - '0') * pow(10, power);
		power--;
	}
	return sum;

}

double sum(char arr[], int start, int end)
{
	double sum = 0;
	int save = 0;
	if (arr[start] == '-')
	{
		if (start == end)
		{
			return -1;
		}
		else
		{
			for (int i = start + 1; i <= end; i++)
			{
				if (arr[i] == '.')
				{
					save = i;
				}
			}
			if (save == 0)
			{
				sum = sum_int(arr, start + 1, end);
			}
			else
				sum = sum_int(arr, start + 1, save - 1) + decimal_sum(arr, save + 1, end);

			return sum * -1;

		}


	}

	else if (arr[start] == '+')
	{
		if (start == end)
		{
			return 1;
		}
		else
		{
			for (int i = start + 1; i <= end; i++)
			{
				if (arr[i] == '.')
				{
					save = i;
				}
			}

			if (save == 0)
			{
				sum = sum_int(arr, start + 1, end);
			}

			else
				sum = sum_int(arr, start + 1, save - 1) + decimal_sum(arr, save + 1, end);

			return sum;


		}

	}
	else
	{


		for (int i = start; i <= end; i++)
		{
			if (arr[i] == '.')
			{
				save = i;
			}
		}
		if (save == 0)
		{
			sum = sum_int(arr, start, end);
		}
		else
			sum = sum_int(arr, start, save - 1) + decimal_sum(arr, save + 1, end);

		return sum;
	}

}

double deffrientiation(double arr[], int size, double x)  //defrientiate our functions only 
{
	double sum = 0;
	double power = size - 1;


	for (int i = 0; i < size; i++)
	{
		sum += power * arr[i] * pow(x, power - 1);
		power--;
		if (power == 1)
		{
			sum += arr[i + 1];
			break;
		}
	}
	return sum;
}

void long_division(double arr[], double arr2[], int size, double root_1)  //specialized for our case onlyy
{
	arr2[0] = arr[0];

	for (int i = 1; i < size - 1; i++)
	{
		arr2[i] = arr[i] + (arr2[i - 1] * root_1);
	}

}

double newton_raphson(double arr[], int degree, double guess, int counter)
{
	double x1 = 0;
	do
	{
		x1 = guess;
		guess = x1 - calculate(arr, degree, x1) / deffrientiation(arr, degree + 1, x1);
		counter++;
		if (counter == 4000)
		{
			guess = newton_raphson(arr, degree, -guess, 4001);
			break;
		}
		if (counter == 8000)
		{
			guess = 12345;
			break;
		}
	} while (fabs(guess - x1) >= 0.01);

	if (guess != 12345)
	{
		cout << "By Newton Raphson Method" << endl << endl;
	}
	return guess;
}

void second_degree(double arr[], int size)
{
	double c = arr[size - 1];
	double b = arr[size - 2];
	double a = arr[size - 3];

	double z = (pow(b, 2) - 4 * a * c);
	if (z > 0)
	{
		cout << "r1=" << (-b + sqrt(z)) / (2 * a) << endl;
		cout << "r2=" << (-b - sqrt(z)) / (2 * a) << endl;
	}
	else if (z == 0)
	{
		cout << "r1=r2=" << -b / (2 * a) << endl;
	}
	else
	{
		cout << "r1=" << -b / (2 * a) << "+" << sqrt(-z) / (2 * a) << "i" << endl;
		cout << "r2=" << -b / (2 * a) << "-" << sqrt(-z) / (2 * a) << "i" << endl;
	}
}

void third_degree(double arr[], int size, double& root_3, double guess)
{
	double d = arr[size - 1];
	double c = arr[size - 2];
	double b = arr[size - 3];
	double a = arr[size - 4];

	if (d == 0)
		root_3 = 0;

	else if ((a + b + c + d) == 0)
		root_3 = 1;

	else if ((b + d) == (a + c))
		root_3 = -1;
	else                                      // searching for real integer roots first 
	{
		int A_counter = 1, D_counter = 1;

		for (int i = 2; i <= a; i++)
		{
			if (int(a) % i == 0)
			{
				A_counter++;
			}
		}
		for (int i = 2; i <= d; i++)
		{
			if (int(d) % i == 0)
			{
				D_counter++;
			}
		}
		D_counter *= 2;
		//(3x ^ 3 + 2x ^ 2 + 2x + 9)

		int* A_arr = new int[A_counter];
		int* D_arr = new int[D_counter];

		D_arr[0] = 1, D_arr[1] = -1;
		A_arr[0] = 1;

		int z = 1;
		for (int i = 2; i <= a; i++)
		{
			if (int(a) % i == 0)
			{
				A_arr[z] = i;
				z++;
			}
		}
		z = 2;
		for (int i = 2; i <= d; i++)
		{
			if (int(d) % i == 0)
			{
				D_arr[z] = i;
				D_arr[z + 1] = -i;
				z += 2;
			}
		}


		double* AD_arr = new double[A_counter * D_counter];

		z = 0;
		for (int k = 0; k < A_counter; k++)
		{
			for (int i = 0; i < D_counter; i++)
			{
				AD_arr[z] = double(D_arr[i]) / double(A_arr[k]);
				z++;
			}
		}



		for (int i = 0; i < A_counter * D_counter; i++)
		{
			if (calculate(arr, 3, AD_arr[i]) == 0)
			{
				root_3 = AD_arr[i];
				break;
			}
		}
	}

	if (root_3 == 12345)
	{
		root_3 = newton_raphson(arr, 3, (d / a), 0);
	}


	if (root_3 == 12345)
	{
		double q, z, t, u1, u2, n;
		q = 2 * pow(b, 3) - 9 * a * b * c + 27 * pow(a, 2) * d;
		z = 4 * pow((pow(b, 2) - 3 * a * c), 3);
		t = sqrt(pow(q, 2) - z);
		u1 = cbrt(0.5 * (q + t));
		u2 = cbrt(0.5 * (q - t));
		n = 1 / (3 * a);
		root_3 = (-b / (3 * a)) - (n * u1) - (n * u2);
		cout << "By the cubic equation rule " << endl << endl;
	}

	if (root_3 == 12345)       //if no real INTEGER roots we search number by number "root_3 will stay the same"
	{
		double save = 0;
		double min = fabs(calculate(arr, 3, -1000000));

		for (int i = -1000000; i < 1000000; i += 500)
		{
			if (i == 0)
				continue;

			if (fabs(calculate(arr, 3, i)) < min)
			{

				min = fabs(calculate(arr, 3, i));
				save = i;
			}
		}

		cout << endl;
		cout << save << endl;

		if (save > 0)
		{
			save += 500;
			double z = (save / 4000);

			for (double i = save; i > 0; i -= z)
			{
				if (fabs(calculate(arr, 3, i)) < min && (i / save) < 1.2)
				{
					min = fabs(calculate(arr, 3, i));
					save = i;
				}
			}
			save += z;

			for (double i = save; i > 0; i -= 0.001)
			{
				if (fabs(calculate(arr, 3, i)) < min)
				{
					min = fabs(calculate(arr, 3, i));
					root_3 = i;
				}
			}
		}
		else if (save < 0)
		{
			save -= 500;
			double z = (-save / 4000);

			for (double i = save; i < 0; i += z)
			{
				if (fabs(calculate(arr, 3, i)) < min && (save / i) < 1.2)
				{
					min = fabs(calculate(arr, 3, i));
					save = i;
					cout << save << endl;
				}
			}
			cout << endl;
			save -= z;
			cout << save << endl;
			for (double i = save; i < 0; i += 0.001)
			{
				if (fabs(calculate(arr, 3, i)) < min)
				{
					min = fabs(calculate(arr, 3, i));
					root_3 = i;
				}
			}
		}
	}



	double arr2[3];
	long_division(arr, arr2, 4, root_3);

	second_degree(arr2, 3);

}

void fourth_degree(double arr[], int size, double& root_4, double& root_3, complex <double> roots[4])
{
	double e = arr[size - 1];
	double d = arr[size - 2];
	double c = arr[size - 3];
	double b = arr[size - 4];
	double a = arr[size - 5];

	if (e == 0)
		root_4 = 0;

	else if ((a + b + c + d + e) == 0)
		root_4 = 1;

	else                             // searching for real integer roots first 
	{
		int A_counter = 1, E_counter = 1;

		for (int i = 2; i <= a; i++)
		{
			if (int(a) % i == 0)
			{
				A_counter++;
			}
		}
		for (int i = 2; i <= e; i++)
		{
			if (int(e) % i == 0)
			{
				E_counter++;
			}
		}
		E_counter *= 2;

		int* A_arr = new int[A_counter];
		int* E_arr = new int[E_counter];

		E_arr[0] = 1, E_arr[1] = -1;
		A_arr[0] = 1;

		int z = 1;
		for (int i = 2; i <= a; i++)
		{
			if (int(a) % i == 0)
			{
				A_arr[z] = i;
				z++;
			}
		}
		z = 2;
		for (int i = 2; i <= e; i++)
		{
			if (int(e) % i == 0)
			{
				E_arr[z] = i;
				E_arr[z + 1] = -i;
				z += 2;
			}
		}
		double* AE_arr = new double[A_counter * E_counter];

		z = 0;
		for (int k = 0; k < A_counter; k++)
		{
			for (int i = 0; i < E_counter; i++)
			{
				AE_arr[z] = double(E_arr[i]) / double(A_arr[k]);
				z++;
			}
		}



		for (int i = 0; i < A_counter * E_counter; i++)
		{
			if (calculate(arr, 4, AE_arr[i]) == 0)
			{
				root_4 = AE_arr[i];
				break;
			}
		}
	}


	if (root_4 == 12345)
	{
		root_4 = newton_raphson(arr, 4, e / a, 0);
	}



	if (root_4 == 12345)
	{
		b /= a;
		c /= a;
		d /= a;
		e /= a;

		complex <double> b2 = b * b;
		complex <double> b3 = b * b2;
		complex <double> b4 = b2 * b2;

		complex <double> alpha = (-3.0 / 8.0) * b2 + c;
		complex <double> beta = b3 / 8.0 - b * c / 2.0 + d;
		complex <double> gamma = (-3.0 / 256.0) * b4 + b2 * c / 16.0 - b * d / 4.0 + e;

		complex <double> alpha2 = alpha * alpha;
		complex <double> t = -b / 4.0;

		if (real(beta) == 0)
		{
			complex <double> rad = sqrt(alpha2 - 4.0 * gamma);
			complex <double> r1 = sqrt((-alpha + rad) / 2.0);
			complex <double> r2 = sqrt((-alpha - rad) / 2.0);

			roots[0] = t + r1;
			roots[1] = t - r1;
			roots[2] = t + r2;
			roots[3] = t - r2;
		}
		else
		{
			complex <double> alpha3 = alpha * alpha2;
			complex <double> P = -(alpha2 / 12.0 + gamma);
			complex <double> Q = -alpha3 / 108.0 + alpha * gamma / 3.0 - beta * beta / 8.0;
			complex <double> R = -Q / 2.0 + sqrt(Q * Q / 4.0 + P * P * P / 27.0);
			complex <double> U = cbrt(R, 0);
			complex <double> y = (-5.0 / 6.0) * alpha + U;
			if (real(U) == 0)
			{
				y -= cbrt(Q, 0);
			}
			else
			{
				y -= P / (3.0 * U);
			}
			complex <double> W = sqrt(alpha + 2.0 * y);

			complex <double> r1 = sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / W));
			complex <double> r2 = sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / W));

			roots[0] = t + (W - r1) / 2.0;
			roots[1] = t + (W + r1) / 2.0;
			roots[2] = t + (-W - r2) / 2.0;
			roots[3] = t + (-W + r2) / 2.0;
		}


	}
	if (root_4 != 12345)
	{
		double arr3[4];

		long_division(arr, arr3, 5, root_4);

		third_degree(arr3, 4, root_3, 2);
	}

}

void fifth_degree(double arr[], int size, double& root_5, double& root_4, double& root_3, complex <double> roots[4])
{
	double f = arr[size - 1];
	double e = arr[size - 2];
	double d = arr[size - 3];
	double c = arr[size - 4];
	double b = arr[size - 5];
	double a = arr[size - 6];

	if (f == 0)
		root_5 = 0;

	else if ((a + b + c + d + e + f) == 0)
		root_5 = 1;


	else                       // searching for real integer roots first 
	{
		int A_counter = 1, F_counter = 1;

		for (int i = 2; i <= a; i++)    //counting number of "possible" roots
		{
			if (int(a) % i == 0)
			{
				A_counter++;
			}
		}
		for (int i = 2; i <= f; i++)
		{
			if (int(f) % i == 0)
			{
				F_counter++;
			}
		}
		F_counter *= 2;

		int* A_arr = new int[A_counter];
		int* F_arr = new int[F_counter];

		F_arr[0] = 1, F_arr[1] = -1;
		A_arr[0] = 1;

		int z = 1;
		for (int i = 2; i <= a; i++)
		{
			if (int(a) % i == 0)
			{
				A_arr[z] = i;
				z++;
			}
		}
		z = 2;
		for (int i = 2; i <= f; i++)
		{
			if (int(f) % i == 0)
			{
				F_arr[z] = i;
				F_arr[z + 1] = -i;
				z += 2;
			}
		}
		double* AF_arr = new double[A_counter * F_counter];

		z = 0;
		for (int k = 0; k < A_counter; k++)
		{
			for (int i = 0; i < F_counter; i++)
			{
				AF_arr[z] = double(F_arr[i]) / double(A_arr[k]);
				z++;
			}
		}



		for (int i = 0; i < A_counter * F_counter; i++)
		{

			if (calculate(arr, 5, AF_arr[i]) == 0)
			{
				root_5 = AF_arr[i];
				cout << "By rational root theorem " << endl << endl;
				break;

			}

		}

	}
	if (root_5 == 12345)
	{
		root_5 = newton_raphson(arr, 5, (f / a), 0);
	}


	if (root_5 == 12345)       //if no real INTEGER roots we search number by number "root_3 will stay the same"
	{
		double save = 0;
		double min = fabs(calculate(arr, 5, -1000000));

		for (int i = -1000000; i < 1000000; i += 500)
		{
			if (i == 0)
				continue;

			if (fabs(calculate(arr, 5, i)) < min)
			{

				min = fabs(calculate(arr, 5, i));
				save = i;
			}
		}


		if (save > 0)
		{
			save += 500;
			double z = (save / 4000);

			for (double i = save; i > 0; i -= z)
			{
				if (fabs(calculate(arr, 5, i)) < min && (i / save) < 1.2)
				{
					min = fabs(calculate(arr, 5, i));
					save = i;
				}
			}
			save += z;

			for (double i = save; i > 0; i -= 0.001)
			{
				if (fabs(calculate(arr, 5, i)) < min)
				{
					min = fabs(calculate(arr, 5, i));
					root_5 = i;
				}
			}
		}
		else if (save < 0)
		{
			save -= 500;
			double z = (-save / 4000);

			for (double i = save; i < 0; i += z)
			{
				if (fabs(calculate(arr, 5, i)) < min && (i / save) < 1.2)
				{
					min = fabs(calculate(arr, 5, i));
					save = i;

				}
			}

			save -= z;

			for (double i = save; i < 0; i += 0.001)
			{
				if (fabs(calculate(arr, 3, i)) < min)
				{
					min = fabs(calculate(arr, 3, i));
					root_5 = i;
				}
			}
			cout << "By iterration method " << endl << endl;
		}
	}


	double arr4[5];

	long_division(arr, arr4, 6, root_5);

	fourth_degree(arr4, 5, root_4, root_3, roots);

}

void first_degree(double arr[], int size)
{
	cout << "root is :" << -arr[size - 1] / arr[size - 2];
}

bool error(char str[], int size)
{
	if (str[size - 1] == '+' || str[size - 1] == '-')
	{
		cout << "syntax error" << endl;
		return 1;
	}
	else
	{
		for (int i = 0; i < size; i++)
		{
			if (str[i] == '^' && str[i - 1] != 'x')
			{
				cout << "syntax error: variable missing" << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '+' && (str[i + 1] == '+'))
			{
				cout << "syntax error: \"++\"" << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '-' && (str[i + 1] == '-'))
			{
				cout << "syntax error: \"--\"" << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '+' && (str[i + 1] == '-'))
			{
				cout << "syntax error: \"+-\"" << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '-' && (str[i + 1] == '+'))
			{
				cout << "syntax error: \"-+\"" << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '^' && (str[i + 1] > '5'))
			{
				cout << "syntax error:missing power, or higher degree max (5) " << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '^' && (str[i + 1] <= '0'))
			{
				cout << "syntax error:missing power, or higher degree max (5) " << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '^' && str[i] == '+')
			{
				cout << "syntax error:missing power" << endl << "input again :" << endl;
				return 1;
			}
			else if (str[i] == '^' && str[i] == '-')
			{
				cout << "syntax error:missing power" << endl << "input again :" << endl;
				return 1;

			}
			else if (str[i] == 'x' && str[i + 1] >= '0' && str[i + 1] <= '9')
			{
				cout << "syntax error:missing power sign" << endl << "input again :" << endl;
				return 1;

			}
		}
		/*
		for (int k=0; k < size; k++)
		{
			if (str[k] == 'x')
				break;
		}
		if (k == size - 1)
		{
			cout << "syntax error:no 'x'" << endl << "input again :" << endl;
			return 1;
		}
	*/
	}
	return 0;
}

int missing_term(char str[], int max)
{
	// 25x ^ 3 - 26x ^ 4 - 23x ^ 2 - 2x + 1 +x^5     max=5     

	int counter = 1;

	for (int i = 2; str[i] != '\0'; i++)
	{
		if (str[i] == '+' || str[i] == '-')
		{
			counter++;
		}
	}
	return counter - max - 1;
}

void call(char str[], double& root_5, double& root_4, double& root_3, int size, double guess, complex <double> roots[4])
{

	int max = 0;  // max will be equal number of terms - 1
	for (int i = 0; i < size; i++)
	{
		if (str[i] == '^' && (str[i + 1] - '0') > max)  // max is equal to the greatest power , number of terms is > than highest power by 1
		{
			max = str[i + 1] - '0';
		}
	}
	if (max == 0)   // if no ^..  the number of terms is 2 ex. "3x+1" no power 
	{
		max++;
	}

	double* arr = new double[max + 1]; // dynamic array of size maximum power + 1
	for (int i = 0; i < max + 1; i++)
	{
		arr[i] = 0;
	}

	int i = 0;

	if (str[0] == '-' || str[0] == '+')
		i++;




	int save = 0;



	for (; str[i] != '+' && str[i] != '-'; i++)   //first term input
	{

		if (str[i] == 'x' && str[i + 1] == '^')
		{
			save = str[i + 2] - '0';
			if (i - 1 != -1)
			{
				arr[max - save] = sum(str, 0, i - 1);
				break;

			}
			else if (i - 1 == -1)
			{
				arr[max - save] = sum(str, 0, 0);
			}

		}
		if (str[0] == 'x')
		{
			arr[max - save] = 1;
		}


		else if (str[i] == 'x' && (str[i + 1] == '+' || str[i + 1] == '-'))
		{
			save = 1;
			arr[max - save] = sum(str, 0, i - 1);
			break;

		}
	}

	if (str[i] == '+' || str[i] == '-')
	{
		arr[max] = sum(str, 0, i - 1);
	}


	int j = i;
	i = size - 1;
	int k = 0;

	for (; str[i] != '+' && str[i] != '-'; i--)   //last term input 
	{
	}

	int m = i;
	i = size - 1;

	for (; str[i] != '+' && str[i] != '-'; i--)   //last term input 
	{

		if (str[i] == 'x' && str[i + 1] == '^')
		{
			save = str[i + 2] - '0';                     //2x + 12x ^ 2 + 30x ^ 3 + 5 + 44x ^ 4
			arr[max - save] = sum(str, m, i - 1);

			break;
		}
		else if (str[i] == 'x' && (i + 1) == size)
		{
			save = 1;
			arr[max - save] = sum(str, m, i - 1);

			break;
		}
	}

	if (str[i] == '+' || str[i] == '-')
	{
		arr[max] = sum(str, i, size - 1);
	}

	int start = 0;
	int end = 0;

	for (; j < m; j++)    //middle terms input              3x^3+25x^2-69x+420
	{
		if (str[j] == '+' || str[j] == '-')
		{
			start = j;
			j++;

			for (; str[j] != '+' && str[j] != '-'; j++)
			{
				if (str[j] == 'x' && str[j + 1] == '^')
				{
					save = str[j + 2] - '0';
					arr[max - save] = sum(str, start, j - 1);
					break;
				}
				else if (str[j] == 'x' && (str[j + 1] == '+' || str[j + 1] == '-'))
				{
					save = 1;
					arr[max - save] = sum(str, start, j - 1);
					break;

				}
			}
			if (str[j] == '+' || str[j] == '-')
			{
				arr[max] = sum(str, start, j - 1);
			}
		}
	}




	cout << endl;
	for (int i = 0; i < max + 1; i++)
	{
		cout << arr[i] << "\t";
	}
	cout << endl << endl;



	switch (max)
	{
	case 5:
		fifth_degree(arr, max + 1, root_5, root_4, root_3, roots);
		break;
	case 4:
		fourth_degree(arr, max + 1, root_4, root_3, roots);
		break;
	case 3:
		third_degree(arr, max + 1, root_3, guess);
		break;
	case 2:
		second_degree(arr, max + 1);
		break;
	case 1:
		first_degree(arr, max + 1);
	default:
		break;
	}
}





int main()
{
	cout << "Enter Polynomial Function: Ax^5+Bx^4+Cx^3+Dx^2+Ex+F\n";
	string x;
	getline(cin, x);
	int space_count = 0;
	for (int i = 0; i < x.length(); i++)
	{
		if (x[i] == ' ')
		{
			space_count++;
		}
	}

	char* str = new char[x.length() + 1 - space_count];

	int size = x.length() - space_count;
	int j = 0;
	for (int i = 0; i < x.length() + 1; i++)
	{
		if (x[i] == ' ')
			continue;

		str[j] = x[i];
		j++;
	}

	int z = error(str, size);


	/*while (z == 1)
	{
		getline(cin, x);
		int space_count = 0;
		for (int i = 0; i < x.length(); i++)
		{
			if (x[i] == ' ')
			{
				space_count++;
			}
		}

		char* str = new char[x.length() + 1 - space_count];

		int size = x.length() - space_count;
		int j = 0;
		for (int i = 0; i < x.length() + 1; i++)
		{
			if (x[i] == ' ')
				continue;

			str[j] = x[i];
			j++;
		}
		z=error(str, size);
	}
	*/




	complex <double> roots[4];

	double root_3 = 12345;  //third degree

	double root_4 = 12345;  //fourth degree

	double root_5 = 12345;  //fifth degree
	if (z == 0)
	{
		call(str, root_5, root_4, root_3, size, 0, roots);

		if (root_3 != 12345)
			cout << "r3 =" << root_3 << endl;

		if (root_4 != 12345)
			cout << "r4 =" << root_4 << endl;

		if (root_5 != 12345)
			cout << "r5 =" << root_5 << endl;
	}
	if (root_3 == 12345)
	{
		for (int i = 0; i < 4; i++)
		{

			if (imag(roots[i]) != 0)

				cout << "R" << i + 1 << ": " << real(roots[i]) << "+" << imag(roots[i]) << "i" << endl;


			else

				cout << "R" << i + 1 << ": " << real(roots[i]) << endl;

		}

	}

	return 0;
}