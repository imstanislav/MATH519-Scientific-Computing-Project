#include <bits/stdc++.h>
#define w1 freopen("v_star.txt","w", stdout)
#define w2 freopen("u_star.txt","w", stdout)
#define w3 freopen("u.txt","w", stdout)
#define w4 freopen("p.txt","w", stdout)
#define w5 freopen("v.txt","w", stdout)
#define w6 freopen("bbfsol.dat","w", stdout)
#define w7 freopen("profile.txt","w", stdout)
#define pi 3.1415926535897932
using namespace std;
const double X = 1.0, Y = 1.0, dt = 1e-5, dx = 0.01, dy = 0.01, visc = 1.0 / 1.0, ro = 1.0, epsP = 1e-7, eps = 1e-3;
const double dxdx = dx * dx, dydy = dy * dy, c1 = ro * dxdx * dydy / (-2 * (dxdx + dydy) * dt), c2 = dydy / (-2 * (dxdx + dydy)), c3 = dxdx / (-2 * (dxdx + dydy));
const double newC = 2 / dxdx + 2 / dydy;
const int x_size = X / dx + 1, y_size = Y / dy + 1, iter_num = 100000;
double u[2][x_size][y_size], v[2][x_size][y_size], u_star[x_size][y_size], v_star[x_size][y_size], P[2][x_size][y_size], maxerrP, temperrP, maxerr, temperr;
bool mode = false;
int curP = 1, prevP = 0;
int cur = 1, prevU = 0;
void showsol() {
	w6;
	cout << "VARIABLES=" << "\"X\"," << "\"Y\"," << "\"U\"," << "\"V\"," << "\"P\"" << endl;
	for (int i = x_size - 1; i >= 0; i--) {
		for (int j = 0; j < y_size; j++) {
			printf("%lf %lf %lf %lf %lf\n", i * dx, j * dy, u[prevU][i][j], v[prevU][i][j], P[prevP][i][j]);
		}
	}
	printf("\n");
}
void showvstar() {
	w1;
	for (int j = x_size - 1; j >= 0; j--) {
		for (int i = 0; i < x_size; i++) {
			printf("%lf \t", v_star[i][j]);
		}
		printf("\n");
	}
}
void showustar() {
	w2;
	for (int j = x_size - 1; j >= 0; j--) {
		for (int i = 0; i < x_size; i++) {
			printf("%lf \t", u_star[i][j]);
		}
		printf("\n");
	}
}
void showp() {
	w4;
	for (int j = x_size - 1; j >= 0; j--) {
		for (int i = 0; i < x_size; i++) {
			printf("%lf \t", P[prevP][i][j]);
		}
		printf("\n");
	}
}
void showu() {
	w3;
	for (int j = x_size - 1; j >= 0; j--) {
		for (int i = 0; i < x_size; i++) {
			printf("%lf \t", u[prevU][i][j]);
		}
		printf("\n");
	}
}
void showv() {
	w5;
	for (int j = x_size - 1; j >= 0; j--) {
		for (int i = 0; i < x_size; i++) {
			printf("%lf \t", v[prevU][i][j]);
		}
		printf("\n");
	}
}

double T(double t) {
	return 1 - exp(-t);
}

double analyticu(double x, double y, double t) {
	return T(t) * sin(pi * x) * cos(pi * y);
}

double analyticv(double x, double y, double t) {
	return - T(t)* cos(pi * x) * sin(pi * y);
}

double analyticp(double x, double y, double t) {
	return 0;
}

double initialu(double x, double y) {
	return 0;
}
double initialv(double x, double y) {
	return 0;
}
double initialp(double x, double y) {
	return 0;
}

//double f1(double x, double y) {
//	return pi * (16 * pi * visc * sin(pi * x) * sin(pi * x) * sin(2 * pi * y) - sin(pi * x) * cos(pi * x) * sin(2 * pi * y) - 2 * sin(pi * y) * sin(pi * y) * cos(2 * pi * x) * sin(pi * x) - sin(pi * y) * cos(pi * x));
//}

//double f1(double x, double y, double t) {
//	return 0.5 * pi * pi * pi * Et(t) * Et(t) * sin(2 * pi * x);
//}

double f1(double x, double y, double t) {
	return sin(pi * x) * cos(pi * y) * (2 * pi * pi + exp(-t) * (1 - 2 * pi * pi)) + 0.5 * pi * sin(2 * pi * x) * (1 - exp(-t)) * (1 - exp(-t));
}

//double f2(double x, double y) {
//	return pi * (-16 * pi * visc * sin(pi * y) * sin(pi * y) * sin(2 * pi * x) + sin(pi * y) * cos(pi * y) * sin(2 * pi * x) + 2 * sin(pi * x) * sin(pi * x) * cos(2 * pi * y) * sin(pi * y) + sin(pi * x) * cos(pi * y));
//}

//double f2(double x, double y, double t) {
//	return 0.5 * pi * pi * pi * Et(t) * Et(t) * sin(2 * pi * y);
//}

double f2(double x, double y, double t) {
	return -cos(pi * x) * sin(pi * y) * (2 * pi * pi + exp(-t) * (1 - 2 * pi * pi)) + 0.5 * pi * sin(2 * pi * y) * (1 - exp(-t)) * (1 - exp(-t));
}

int main() {
	memset(u, sizeof(u), 0);
	memset(v, sizeof(v), 0);
	memset(u_star, sizeof(u_star), 0);
	memset(v_star, sizeof(v_star), 0);
	memset(P, sizeof(P), 0);
	//boundaries and inits.
	//inits
	for (int i = 0; i < x_size; i++) {
		for (int j = 0; j < y_size; j++) {
			u[cur][i][j] = initialu(i * dx, j * dy);
			u[prevU][i][j] = initialu(i * dx, j * dy);

			v[cur][i][j] = initialv(i * dx, j * dy);
			v[prevU][i][j] = initialv(i * dx, j * dy);

			P[curP][i][j] = initialp(i * dx, j * dy);
			P[prevP][i][j] = initialp(i * dx, j * dy);

		}
	}

	//boundary
	int i, j, iter = 0;

	do {

		//update boundaries
		for (int i = 0; i < x_size; i++) {
			u[cur][i][0] = analyticu(i * dx, 0, iter * dt);
			u[prevU][i][0] = analyticu(i * dx, 0, iter * dt);
			u_star[i][0] = analyticu(i * dx, 0, iter * dt);

			u[cur][i][y_size - 1] = analyticu(i * dx, 1, iter * dt);
			u[prevU][i][y_size - 1] = analyticu(i * dx, 1, iter * dt);
			u_star[i][y_size - 1] = analyticu(i * dx, 1, iter * dt);


			v[cur][i][0] = analyticv(i * dx, 0, iter * dt);
			v[prevU][i][0] = analyticv(i * dx, 0, iter * dt);
			v_star[i][0] = analyticv(i * dx, 0, iter * dt);


			v[cur][i][y_size - 1] = analyticv(i * dx, 1, iter * dt);
			v[prevU][i][y_size - 1] = analyticv(i * dx, 1, iter * dt);
			v_star[i][y_size - 1] = analyticv(i * dx, 1, iter * dt);
		}


		for (int j = 0; j < y_size; j++) {
			u[cur][0][j] = analyticu(0, j * dy, iter * dt);
			u[prevU][0][j] = analyticu(0, j * dy, iter * dt);
			u_star[0][j] = analyticu(0, j * dy, iter * dt);

			u[cur][x_size - 1][j] = analyticu(1, j * dy, iter * dt);
			u[prevU][x_size - 1][j] = analyticu(1, j * dy, iter * dt);
			u_star[x_size - 1][j] = analyticu(1, j * dy, iter * dt);


			v[cur][0][j] = analyticv(0, j * dy, iter * dt);
			v[prevU][0][j] = analyticv(0, j * dy, iter * dt);
			v_star[0][j] = analyticv(0, j * dy, iter * dt);


			v[cur][x_size - 1][j] = analyticv(1, j * dy, iter * dt);
			v[prevU][x_size - 1][j] = analyticv(1, j * dy, iter * dt);
			v_star[x_size - 1][j] = analyticv(1, j * dy, iter * dt);
		}

		/*u and v stars*/
		temperr = -99999;
		for (i = 1; i < x_size - 1; i++) {
			for (j = 1; j < y_size - 1; j++) {
				u_star[i][j] = u[prevU][i][j] + dt * (f1(i*dx, j*dy, iter * dt) + visc * ((u[prevU][i + 1][j] - 2 * u[prevU][i][j] + u[prevU][i - 1][j]) / dxdx + (u[prevU][i][j + 1] - 2 * u[prevU][i][j] + u[prevU][i][j - 1]) / dydy)
					- u[prevU][i][j] * (u[prevU][i][j] - u[prevU][i - 1][j]) / (dx)-v[prevU][i][j] * (u[prevU][i][j] - u[prevU][i][j - 1]) / (dy));
				v_star[i][j] = v[prevU][i][j] + dt * (f2(i *dx, j *dy, iter * dt) + visc * ((v[prevU][i + 1][j] - 2 * v[prevU][i][j] + v[prevU][i - 1][j]) / dxdx + (v[prevU][i][j + 1] - 2 * v[prevU][i][j] + v[prevU][i][j - 1]) / dydy)
					- u[prevU][i][j] * (v[prevU][i][j] - v[prevU][i - 1][j]) / (dx)-v[prevU][i][j] * (v[prevU][i][j] - v[prevU][i][j - 1]) / (dy));
			}
		}
		//neuman upper wall
		//for (i = 0; i < x_size; i++) {
		//	v_star[i][y_size - 1] = v_star[i][y_size - 2];
		//}
		/*pressure*/
		do {
			temperrP = -9999999;
			for (i = 1; i < x_size - 1; i++) {
				for (j = 1; j < y_size - 1; j++) {
					P[curP][i][j] = -dt / ro * ((u_star[i][j] - u_star[i - 1][j]) / (dx)+(v_star[i][j] - v_star[i][j - 1]) / (dy))/newC
						- (P[prevP][i + 1][j] + P[prevP][i - 1][j])/(dxdx*newC) - (P[prevP][i][j + 1] + P[prevP][i][j - 1])/(dydy * newC);
					temperrP = max(temperrP, abs(P[curP][i][j] - P[prevP][i][j]));
				}
			}
			//neuman
			maxerrP = temperrP;
			curP = !curP;
			prevP = !prevP;
			//printf("%lf\n", maxerrP);
		} while (maxerrP > epsP);
		for (i = 1; i < x_size - 1; i++) {
			for (j = 1; j < y_size - 1; j++) {
				u[cur][i][j] = -dt / ro * (P[prevP][i + 1][j] - P[prevP][i][j]) / (dx)+u_star[i][j];
				v[cur][i][j] = -dt / ro * (P[prevP][i][j + 1] - P[prevP][i][j]) / (dy)+v_star[i][j];
				temperr = max(temperr, abs(u[cur][i][j] - analyticu(i*dx, j*dy, iter * dt)));
				temperr = max(temperr, abs(v[cur][i][j] - analyticv(i * dx, j * dy, iter * dt)));

			}
		}
		curP = !curP;
		cur = !cur;
		prevP = !prevP;
		prevU = !prevU;
		maxerr = temperr;
		//neuman upper wall
		//for (i = 0; i < x_size; i++) {
		//	v[i][y_size - 1] = v[i][y_size - 2];
		//}
		printf("%d %lf\n", iter, maxerr);
		iter++;
	} while (iter < iter_num);

	printf("Test 1: Maximal error is %lf after %d iterations", maxerr, iter_num);
	showp();
	showu();
	showv();
	showustar();
	showvstar();
	showsol();
	//w7;
	//for (int j = 0; j < y_size; j++) {
	//	printf("%lf %lf\n", j * dy, u[24][j]);
	//}
	//return 0;
}
