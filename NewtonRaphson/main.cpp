#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <functional>
using namespace std;

const double PI = 3.141592653589793;
const int Num = 4;

class NewtonRaphson
{
public:
	function<double(double)> func;
	double X_0 = 0.;
	double X_1 = 0.;
	double dx	= 0.0001;
	double eta = 0.0001;
	double EE = 100.;
	int iteration = 0;

	NewtonRaphson(function<double(double)> F, double x) : func(F), X_0(x) {}

	double dist() { return abs(X_0 - X_1); }
	double dfdx(double x, double dx) { return (func(x + dx) - func(x)) / dx; }
	void excute()
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			X_1 = X_0;
			X_0 = perform();
			EE = dist();
		}
		cout << "Iteration : " << iteration << endl;
		cout << "Answer is : " << X_0 << " And EE : " << EE << endl;
	}
private:
	double perform() { return X_0 - func(X_0) / dfdx(X_0, dx); }
};

int main(int argc, char* argv[])
{
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> dist(0, 10);

	function<double(double)> func([](double x)->double { return (x + 2) * (x - 4); });

	NewtonRaphson NR(func, dist(gen));
	cout << "Starting At : " << NR.X_0 << endl;
	NR.excute();

	cout << endl;

	NewtonRaphson NR2(func, -dist(gen));
	cout << "Starting At : " << NR2.X_0 << endl;
	NR2.excute();

	system("pause");
	return 0;
}