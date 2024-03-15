#pragma once
#include "Data.h"

class UnsteadyFiniteDifferenceExample {
	double rho;
	double u;
	double Gamma;
	int N;
	double dx;
	double dt;
	double d;
	double c;
	double tol;
	std::string temporal;
	void initialize(Data& data) {
		data.x.resize(N + 1);
		for (int i = 0; i < N + 1; i++) {
			data.x[i] = i * dx;
		}
		data.t.push_back(0.0);
		data.phi.resize(1);
		data.phi[0].resize(N + 1);
	}
	void discretize(Data& data, std::vector<double>& A_W, std::vector<double>& A_P, std::vector<double>& A_E, std::vector<double>& Q_P, bool flag) {
		if (temporal == "explicit Euler") {
			for (int i = 0; i < N - 1; i++) {
				A_P[i] = 1;
				Q_P[i] = (d + c / 2) * data.phi[data.phi.size() - 1][i] + (1 - 2 * d) * data.phi[data.phi.size() - 1][i + 1] + (d - c / 2) * data.phi[data.phi.size() - 1][i + 2];
			}
		}
		else if (temporal == "implicit Euler" || (temporal == "three-time-level" && flag)) {
			for (int i = 0; i < N - 1; i++) {
				A_W[i] = -c / 2 - d;
				A_P[i] = 1 + 2 * d;
				A_E[i] = c / 2 - d;
				Q_P[i] = data.phi[data.phi.size() - 1][i + 1];
			}
		}
		else if (temporal == "Crank-Nicolson") {
			for (int i = 0; i < N - 1; i++) {
				A_W[i] = -c / 4 - d / 2;
				A_P[i] = 1 + d;
				A_E[i] = c / 4 - d / 2;
				Q_P[i] = (d / 2 + c / 4) * data.phi[data.phi.size() - 1][i] + (1 - d) * data.phi[data.phi.size() - 1][i + 1] + (d / 2 - c / 4) * data.phi[data.phi.size() - 1][i + 2];
			}
		}
		else if (temporal == "three-time-level" && !flag) {
			for (int i = 0; i < N - 1; i++) {
				A_W[i] = -c / 2 - d;
				A_P[i] = 1.5 + 2 * d;
				A_E[i] = c / 2 - d;
				Q_P[i] = 2 * data.phi[data.phi.size() - 1][i + 1] - 0.5 * data.phi[data.phi.size() - 2][i + 1];
			}
		}
		Q_P[N - 2] -= A_E[N - 2];
	}
	void Thomas(std::vector<double>& L, std::vector<double>& D, std::vector<double>& U, std::vector<double>& b) {
		for (int i = 0; i < D.size() - 1; i++) {
			U[i] = U[i] / D[i];
			D[i + 1] = D[i + 1] - L[i + 1] * U[i];
		}
		b[0] = b[0] / D[0];
		for (int i = 1; i < D.size(); i++) {
			b[i] = (b[i] - L[i] * b[i - 1]) / D[i];
		}
		for (int i = D.size() - 2; 0 <= i; i--) {
			b[i] = b[i] - U[i] * b[i + 1];
		}
	}
	void save(Data& data, std::vector<double>& Q_P) {
		data.t.push_back(data.t[data.t.size() - 1] + dt);
		data.phi.resize(data.phi.size() + 1);
		data.phi[data.phi.size() - 1].resize(N + 1);
		for (int i = 1; i < N; i++) {
			data.phi[data.phi.size() - 1][i] = Q_P[i - 1];
		}
		data.phi[data.phi.size() - 1][N] = 1;
	}
	double residual(Data& data) {
		double res = 0.0;
		for (int i = 0; i < N + 1; i++) {
			res = res < abs(data.phi[data.phi.size() - 1][i] - data.phi[data.phi.size() - 2][i]) ? abs(data.phi[data.phi.size() - 1][i] - data.phi[data.phi.size() - 2][i]) : res;
		}
		return res;
	}
public:
	UnsteadyFiniteDifferenceExample(double rho, double u, double Gamma, int N, double dt, double tol, std::string temporal) :rho(rho), u(u), Gamma(Gamma), N(N), dx(1.0 / N), dt(dt), d(Gamma * dt / (rho * pow(dx, 2))), c(u * dt / dx), tol(tol), temporal(temporal) {}
	void solve(Data& data) {
		std::vector<double> A_W(N - 1), A_P(N - 1), A_E(N - 1), Q_P(N - 1);
		initialize(data);
		bool flag = true;
		while (true) {
			discretize(data, A_W, A_P, A_E, Q_P, flag);
			flag = false;
			Thomas(A_W, A_P, A_E, Q_P);
			save(data, Q_P);
			double res = residual(data);
			if (res < tol) {
				break;
			}
		}
	}
};
