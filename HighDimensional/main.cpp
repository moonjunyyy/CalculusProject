#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

const double PI = 3.141592653589793;

double fxyz(double x, double y, double z)
{
	return sin(PI * x) * sin(PI * y) * sin(PI * z);
}
double dfxyzdx(double x, double y, double z, double dx)
{
	return (fxyz(x + dx, y, z) - fxyz(x, y, z)) / dx;
}
double dfxyzdy(double x, double y, double z, double dy)
{
	return (fxyz(x, y + dy, z) - fxyz(x, y, z)) / dy;
}
double dfxyzdz(double x, double y, double z, double dz)
{
	return (fxyz(x, y, z + dz) - fxyz(x, y, z)) / dz;
}
double distance(double x, double y, double z, double x_1, double y_1,double z_1)
{
	return sqrt((x - x_1) * (x - x_1) + (y - y_1) * (y - y_1) + (z - z_1) * (z - z_1));
}

class GradientDesend3D
{
public:
	double psi = 0.01, eta = 0.0001;
	int iteration = 0;
	double X_0 = 0, Y_0 = 0, Z_0 = 0., X_1 = 0, Y_1 = 0, Z_1 = 0.;
	double EE = 100.;
	double dx = 0.001, dy = 0.001, dz = 0.001;

	GradientDesend3D() { X_0 = 0; Y_0 = 0, Z_0 = 0; }
	GradientDesend3D(double X, double Y, double Z) : X_0(X), Y_0(Y), Z_0(Z) {}

	void setInitCoord(double X, double Y, double Z) { this->X_0 = X, this->Y_0 = Y; this->Z_0 = Z; }
	void initGradient(double X, double Y, double Z)
	{
		setInitCoord(X, Y, Z);
		EE = 100., iteration = 0, X_1 = 0, Y_1 = 0;
	}
	void findLocalMaxima(bool COUT)
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			if (COUT) cout << "f(" << setprecision(4) << X_0 << ", " << setprecision(4) << Y_0 << ", " << setprecision(4) << Z_0 << ") = " << setprecision(6) << fxyz(X_0, Y_0, Z_0) << endl;
			EE = distance(X_0, Y_0, Z_0, X_1, Y_1, Z_1);
			if (COUT) cout << "EE : " << EE << "  diff : " << dfxyzdx(X_0, Y_0, Z_0, dx) << ", " << dfxyzdy(X_0, Y_0, Z_0, dy) << ", " << dfxyzdz(X_0, Y_0, Z_0, dz) << endl;
			X_1 = X_0;
			Y_1 = Y_0;
			Z_1 = Z_0;
			X_0 += psi * dfxyzdx(X_0, Y_0, Z_0, dx);
			Y_0 += psi * dfxyzdy(X_0, Y_0, Z_0, dy);
			Z_0 += psi * dfxyzdz(X_0, Y_0, Z_0, dz);
		}
		cout << endl << iteration << "-th E = " << EE << " f(" << setprecision(6) << X_0 << ", " << setprecision(6) << Y_0 << ", " << setprecision(6) << Z_0 << ")" << endl << endl;
	}
	void findLocalMinima(bool COUT)
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			if (COUT) cout << "f(" << setprecision(4) << X_0 << ", " << setprecision(4) << Y_0 << ", " << setprecision(4) << Z_0 << ") = " << setprecision(6) << fxyz(X_0, Y_0, Z_0) << endl;
			EE = distance(X_0, Y_0, Z_0, X_1, Y_1, Z_1);
			if (COUT) cout << "EE : " << EE << "  diff : " << dfxyzdx(X_0, Y_0, Z_0, dx) << ", " << dfxyzdy(X_0, Y_0, Z_0, dy) << ", " << dfxyzdz(X_0, Y_0, Z_0, dz) << endl;
			X_1 = X_0;
			Y_1 = Y_0;
			Z_1 = Z_0;
			X_0 -= psi * dfxyzdx(X_0, Y_0, Z_0, dx);
			Y_0 -= psi * dfxyzdy(X_0, Y_0, Z_0, dy);
			Z_0 -= psi * dfxyzdz(X_0, Y_0, Z_0, dz);
		}
		cout << endl << iteration << "-th E = " << EE << " f(" << setprecision(6) << X_0 << ", " << setprecision(6) << Y_0 << ", " << setprecision(6) << Z_0 << ")" << endl << endl;
	}
};

int main(int argc, char* argv[])
{
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> dist(0., 1.);

	double X = dist(gen), Y = dist(gen), Z = dist(gen);
	GradientDesend3D GD3D(X, Y, Z);
	cout << "Starting at ( " << setprecision(4) << GD3D.X_0 << " , " << setprecision(4) << GD3D.Y_0 << " , " << setprecision(4) << GD3D.Z_0 << " )" << endl << endl;
	GD3D.findLocalMaxima(true);
}