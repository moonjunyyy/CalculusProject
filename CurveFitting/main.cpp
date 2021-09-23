#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
using namespace std;

const double PI = 3.141592653589793;
const int Num = 4;

class Point
{
public:
	double X, Y;
	Point() : X(0.), Y(0.) {}
	Point(double x, double y) : X(x), Y(y) {}
	double dist(Point P) { return sqrt((X - P.X) * (X - P.X) + (Y - P.Y) * (Y - P.Y)); }
	double dist(double x, double y) { return sqrt((X - x) * (X - x) + (Y - y) * (Y - y)); }
	friend ostream& operator<<(ostream& os, const Point& p)
	{
		os << "(" << p.X << ", " << p.Y << ")";
		return os;
	}
};
class Straight
{
public:
	double A, B, C;
	Straight(double a, double b, double c) : A(a), B(b), C(c) {}
	Point footofPerpendicular(Point P)
	{
		return Point(B * (B * P.X - A * P.Y) - A * C / (A * A + B * B), A * (-B * P.X + A * P.Y) - B * C / (A * A + B * B));
	}
	Point footofPerpendicular(double x, double y)
	{
		return Point(B * (B * x - A * y) - A * C / (A * A + B * B), A * (-B * x + A * y) - B * C / (A * A + B * B));
	}
	friend ostream& operator<<(ostream& os, const Straight& s)
	{
		os << s.A << "x + " << s.B << "y + " << s.C << " = 0";
		return os;
	}
	double sigmaDsquare(float* x, float* y, int N)
	{
		double Sigma = 0;
		for (int i = 0; i < N; i++)
		{
			double D = this->footofPerpendicular(x[i], y[i]).dist(x[i], y[i]);
			Sigma += pow(D,2.0);
		}
		return Sigma;
	}
	double sigmaDsquare(Point* p, int N) 
	{
		double Sigma = 0;
		for (int i = 0; i < N; i++)
		{
			double D = this->footofPerpendicular(p[i]).dist(p[i]);
			Sigma += pow(D, 2.0);
		}
		return Sigma;
	}
};

double dfabda(double a, double b, double da, Point* p, int N)
{
	Straight S1(-a, 1., -b), S2(-a - da, 1., -b);
	return (S2.sigmaDsquare(p, N) - S1.sigmaDsquare(p, N)) / da;
}
double dfabdb(double a, double b, double db, Point* p, int N)
{
	Straight S1(-a, 1., -b), S2(-a, 1., -b - db);
	return (S2.sigmaDsquare(p, N) - S1.sigmaDsquare(p, N)) / db;
}

class GradientDesend2D
{
public:
	double psi = 0.0005, eta = 0.000001;
	int iteration = 0;
	double A_0 = 0, B_0 = 0, A_1 = 0, B_1 = 0;
	double EE = 100.;
	double da = 0.001, db = 0.001;
	Point* points;

	GradientDesend2D() { A_0 = 0; B_0 = 0; points = nullptr; };
	GradientDesend2D(double A, double B, Point* p) :A_0(A), B_0(B), points(p) {}

	void setInitCoord(double A, double B) { this->A_0 = A, this->B_0 = B; }
	void initGradient(double A, double B)
	{
		setInitCoord(A, B);
		EE = 100., iteration = 0, A_1 = 0, B_1 = 0;
	}
	double distance(double x, double y, double x_1, double y_1)
	{
		return sqrt((x - x_1) * (x - x_1) + (y - y_1) * (y - y_1));
	}
	void findLocalMaxima(bool COUT)
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			if (COUT) cout << "y = " << setprecision(4) << A_0 << "x + " << setprecision(4) << B_0 << endl;
			EE = distance(A_0, B_0, A_1, B_1);
			if (COUT) cout << "EE : " << EE << "  diff : " << dfabda(A_0, B_0, da, points, Num) << ", " << dfabdb(A_0, B_0, db, points, Num) << endl;
			A_1 = A_0;
			B_1 = B_0;
			A_0 += psi * dfabda(A_0, B_0, da, points, Num);
			B_0 += psi * dfabdb(A_0, B_0, db, points, Num);
		}
		cout << endl << iteration << "-th E = " << EE << " y = " << setprecision(4) << A_0 << "x + " << setprecision(4) << B_0 << endl;
	}
	void findLocalMinima(bool COUT)
	{
		for (; EE > eta && iteration++ < 10000;)
		{
			if (COUT) cout << "y = " << setprecision(4) << A_0 << "x + " << setprecision(4) << B_0 << endl;
			EE = distance(A_0, B_0, A_1, B_1);
			if (COUT) cout << "EE : " << EE << "  diff : " << dfabda(A_0, B_0, da, points, Num) << ", " << dfabdb(A_0, B_0, db, points, Num) << endl;
			A_1 = A_0;
			B_1 = B_0;
			A_0 -= psi * dfabda(A_0, B_0, da, points, Num);
			B_0 -= psi * dfabdb(A_0, B_0, db, points, Num);
		}
		cout << endl << iteration << "-th E = " << EE << " y = " << setprecision(4) << A_0 << "x + " << setprecision(4) << B_0 << endl;
		
	}
};

int main(int argc, char* argv[])
{
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> dist(-10, 10);
	uniform_real_distribution<double> dist2(-1, 1);

	Point P1(3., 1.);
	Straight S1(1., 0, -1.);

	cout << P1 << endl;
	cout << S1 << endl;

	Point P2 = S1.footofPerpendicular(P1);

	cout << P2 << endl;
	cout << P1.dist(P2) << endl;

	Point* Ps;
	Ps = new Point[Num];
	
	for (int i = 0; i < Num; i++)
	{
		Ps[i].X = dist(gen);
		Ps[i].Y = dist(gen);
	}
	for (int i = 0; i < Num; i++)
	{
		cout << *(Ps + i) << endl;
	}

	GradientDesend2D GD2D(dist2(gen), dist2(gen), Ps);
	for (int i = 0; i < 5; i++)
	{
		cout << "starting at A = " << setprecision(4) << GD2D.A_0 << ", B = " << setprecision(4) << GD2D.B_0;
		GD2D.findLocalMinima(false);
		GD2D.initGradient(dist2(gen), dist2(gen));
		cout << endl;
	}

	double _x = dist2(gen), _y = dist2(gen);
	GradientDesend2D* pGD2D;

	pGD2D = new GradientDesend2D[20];
	for (int i = 0; i < 20; i++)
	{
		pGD2D[i] = GradientDesend2D(_x, _y, Ps);
	}
	for (int i = 0; i < 20; i++)
	{
		pGD2D[i].psi = 0.0005 * (i + 1);
		pGD2D[i].findLocalMinima(false);
	}
	delete[] pGD2D;

	pGD2D = new GradientDesend2D[20];
	for (int i = 0; i < 20; i++)
	{
		pGD2D[i] = GradientDesend2D(_x, _y, Ps);
	}
	for (int i = 0; i < 20; i++)
	{
		pGD2D[i].eta = 0.00001 * (i + 1);
		pGD2D[i].findLocalMinima(false);
	}
	delete[] pGD2D;

	delete[] Ps;

	system("pause");
	return 0;
}