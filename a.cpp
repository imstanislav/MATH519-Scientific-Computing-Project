#include <bits/stdc++.h>
#define w1 freopen("v_star.txt","w", stdout)
#define w2 freopen("u_star.txt","w", stdout)
#define w3 freopen("u.txt","w", stdout)
#define w4 freopen("p.txt","w", stdout)
#define w5 freopen("v.txt","w", stdout)
#define w6 freopen("bbfsol.dat","w", stdout)
#define w7 freopen("profile.txt","w", stdout)
using namespace std;
const double X = 1.0, Y = 1.0, dt = 1e-3, dx = 2 * 1e-2, dy = 2 * 1e-2, visc = 1.0 / 100.0, ro = 1.0, epsP = 1e-4;
const double dxdx = dx * dx, dydy = dy * dy, c1 = ro * dxdx * dydy / (-2 * (dxdx + dydy) * dt), c2 = dydy / (-2 * (dxdx + dydy)), c3 = dxdx / (-2 * (dxdx + dydy));
const int x_size = X / dx + 1, y_size = Y / dy + 1, iter_num = 50000;
double u[x_size][y_size], v[x_size][y_size], u_star[x_size][y_size], v_star[x_size][y_size], P[2][x_size][y_size], maxerrP, temperrP;
bool mode = false;
int curP = 1, prevP = 0;
void showsol() {
	w6;
	cout << "VARIABLES=" << "\"X\"," << "\"Y\"," << "\"U\"," << "\"V\"," << "\"P\"" << endl;
	for (int i = x_size - 1; i >= 0; i--) {
		for (int j = 0; j < y_size; j++) {
			printf("%lf %lf %lf %lf %lf\n", i * dx, j * dy, u[i][j], v[i][j], P[prevP][i][j]);
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
			printf("%lf \t", u[i][j]);
		}
		printf("\n");
	}
}
void showv() {
	w5;
	for (int j = x_size - 1; j >= 0; j--) {
		for (int i = 0; i < x_size; i++) {
			printf("%lf \t", v[i][j]);
		}
		printf("\n");
	}
}
void wantcontinue() {
	double a[5]; char s[20];
	int x, y;
	freopen("bbfsol.dat", "r", stdin);
	scanf("%s", &s);
	for (int i = 0; i < 2601; i++) {
		scanf("%lf %lf %lf %lf %lf\n", &a[0], &a[1], &a[2], &a[3], &a[4]);
		x = a[0] / dx;
		y = a[1] / dy;

		u[x][y] = a[2];
		v[x][y] = a[3];
		P[0][x][y] = a[4];
	}

}
int main() {
	memset(u, sizeof(u), 0);
	memset(v, sizeof(v), 0);
	memset(u_star, sizeof(u_star), 0);
	memset(v_star, sizeof(v_star), 0);
	memset(P, sizeof(P), 0);
	//boundaries and inits.
	for (int i = 0; i < x_size; i++) { // upper wall
		u[i][y_size - 1] = 1.0;
		u_star[i][y_size - 1] = 1.0;
	}
	int i, j, iter = 0;
	if (mode) {
		wantcontinue();
		//w7;
		//for (int j = 0; j < y_size; j++) {
		//	printf("%lf %lf\n", j * dy, u[24][j]);
		//}
		//return 0;
	}

	do {
		/*u and v stars*/
		for (i = 1; i < x_size - 1; i++) {
			for (j = 1; j < y_size - 1; j++) {
				u_star[i][j] = u[i][j] + dt * (visc * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / dxdx + (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / dydy)
					- u[i][j] * (u[i][j] - u[i - 1][j]) / (dx)-v[i][j] * (u[i][j] - u[i][j - 1]) / (dy));
				v_star[i][j] = v[i][j] + dt * (visc * ((v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / dxdx + (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / dydy)
					- u[i][j] * (v[i][j] - v[i - 1][j]) / (dx)-v[i][j] * (v[i][j] - v[i][j - 1]) / (dy));
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
					P[curP][i][j] = c1 * ((u_star[i][j] - u_star[i - 1][j]) / (dx)+(v_star[i][j] - v_star[i][j - 1]) / (dy))
						+ c2 * (P[prevP][i + 1][j] - P[prevP][i - 1][j]) + c3 * (P[prevP][i][j + 1] - P[prevP][i][j - 1]);
					temperrP = max(temperrP, abs(P[curP][i][j] - P[prevP][i][j]));
				}
			}
			//neuman
			for (i = 0; i < x_size; i++) {
				P[curP][i][0] = P[curP][i][1]; // lower wall
				P[curP][i][y_size - 1] = P[curP][i][y_size - 2]; // upper wall
			}
			for (j = 0; j < y_size; j++) {
				P[curP][0][j] = P[curP][1][j]; //left wall
				P[curP][x_size - 1][j] = P[curP][x_size - 2][j]; //right wall
			}
			maxerrP = temperrP;
			curP = !curP;
			prevP = !prevP;
			//printf("%lf\n", maxerrP);
		} while (maxerrP > epsP);
		for (i = 1; i < x_size - 1; i++) {
			for (j = 1; j < y_size - 1; j++) {
				u[i][j] = -dt / ro * (P[prevP][i + 1][j] - P[prevP][i][j]) / (dx)+u_star[i][j];
				v[i][j] = -dt / ro * (P[prevP][i][j + 1] - P[prevP][i][j]) / (dy)+v_star[i][j];
			}
		}
		//neuman upper wall
		//for (i = 0; i < x_size; i++) {
		//	v[i][y_size - 1] = v[i][y_size - 2];
		//}
		printf("%d\n", iter);
		iter++;
	} while (iter < iter_num);
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
