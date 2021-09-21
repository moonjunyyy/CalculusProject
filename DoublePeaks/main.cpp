#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
using namespace std;
const double PI = 3.141592653589793;

//Sprint 1
double gaussian(double x, double y, double mux, double muy, double sigx, double sigy, double peak)
{
	return(peak * exp(-pow((x - mux) / sigx, 2.0) - pow((y - muy) / sigy, 2.0)));
}
double fxy(double x, double y)
{
	return (gaussian(x, y, 1., 1., 1., 2., 4.) + gaussian(x, y, -1., -1., 1., 1., 2.));
}

//Sprint 3
double dfxydx(double x, double y, double dx)
{
	return (fxy(x + dx, y) - fxy(x, y)) / dx;
}
double dfxydy(double x, double y, double dy)
{
	return (fxy(x, y + dy) - fxy(x, y)) / dy;
}
double distance(double x, double y, double x_1, double y_1)
{
	return sqrt((x - x_1) * (x - x_1) + (y - y_1) * (y - y_1));
}

class GradientDesend2D
{
public:
	double psi = 0.01, eta = 0.0001;
	int iteration = 0;
	double X_0 = 0, Y_0 = 0, X_1 = 0, Y_1 = 0;
	double EE = 100.;
	double dx = 0.001, dy = 0.001;

	GradientDesend2D() { X_0 = 0; Y_0 = 0; };
	GradientDesend2D(double X, double Y) :X_0(X), Y_0(Y) {}

	void setInitCoord(double X, double Y) { this->X_0 = X, this->Y_0 = Y; }
	void initGradient(double X, double Y)
	{
		setInitCoord(X, Y);
		EE = 100., iteration = 0, X_1 = 0, Y_1 = 0;
	}
	void findLocalMaxima(bool COUT)
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			if (COUT) cout << "f(" << setprecision(4) << X_0 << ", " << setprecision(4) << Y_0 << ") = " << setprecision(6) << fxy(X_0, Y_0) << endl;
			EE = distance(X_0, Y_0, X_1, Y_1);
			if (COUT) cout << "EE : " << EE << "  diff : " << dfxydx(X_0, Y_0, dx) << ", " << dfxydx(X_0, Y_0, dy) << endl;
			X_1 = X_0;
			Y_1 = Y_0;
			X_0 += psi * dfxydx(X_0, Y_0, dx);
			Y_0 += psi * dfxydy(X_0, Y_0, dy);
		}
		cout << endl << iteration << "-th E = " << EE << "f(" << setprecision(6) << X_0 << ", " << setprecision(6) << Y_0 << ")" << endl << endl;
	}
	void findLocalMinima(bool COUT)
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			if (COUT) cout << "f(" << setprecision(4) << X_0 << ", " << setprecision(4) << Y_0 << ") = " << setprecision(6) << fxy(X_0, Y_0) << endl;
			EE = distance(X_0, Y_0, X_1, Y_1);
			if (COUT) cout << "EE : " << EE << "  diff : " << dfxydx(X_0, Y_0, dx) << ", " << dfxydx(X_0, Y_0, dy) << endl;
			X_1 = X_0;
			Y_1 = Y_0;
			X_0 -= psi * dfxydx(X_0, Y_0, dx);
			Y_0 -= psi * dfxydy(X_0, Y_0, dy);
		}
		cout << endl << iteration << "-th E = " << EE << " at (" << setprecision(6) << X_0 << ", " << setprecision(6) << Y_0 << ")" << endl << endl;
	}
};

int main(int argc, char* argv[])
{
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> dist(-3, 3);

	double X = dist(gen), Y = dist(gen);
	cout << "Starting Point : (" << setprecision(4) << X << ", " << Y << ")" << endl << endl;
	
	GradientDesend2D* GDList;
	GDList = new GradientDesend2D[20];
	cout << endl << endl << "Differnt Psi from 0.001 to 0.4" << endl << endl;

	for (int i = 0; i < 20; i++)
	{
		*(GDList + i) = GradientDesend2D(X, Y);
		(GDList + i)->psi = 0.001 * i * i;
	}
	for (int i = 0; i < 20; i++)
	{
		GDList[i].findLocalMaxima(false);
	}

	delete[] GDList;

	GDList = new GradientDesend2D[20];
	cout << endl << endl << "Differnt eta from 0.0001 to 0.002" << endl << endl;

	for (int i = 0; i < 20; i++)
	{
		*(GDList + i) = GradientDesend2D(X, Y);
		(GDList + i)->eta = 0.0001 * i;
	}
	for (int i = 0; i < 20; i++)
	{
		GDList[i].findLocalMaxima(false);
	}

	delete[] GDList;

	system("pause");
	return 0;
}