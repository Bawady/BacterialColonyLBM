///----------------------------------------------------------------------------------------------------------------------------------
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
// https://www.sciencedirect.com/science/article/pii/S0378437198003458
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot = true;
const double cs2 = 1. / 3.;
const int np = 5;
const vector<int> cx = {0, 1, 0, -1, 0}, cy = {0, 0, 1, 0, -1},
                  opp = {0, 3, 4, 1, 2};
const vector<double> wf = {1 / 3., 1 / 6., 1 / 6., 1 / 6., 1 / 6.};
/*
const int np = 9;
const vector<int> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1},
                                                                  cy = {0, 0, 1,
0, -1, 1, 1, -1, -1}; const vector<double> wf =
{4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
*/
const double L_phys = 1, beta_phys = 0., dN_phys = 1., dB_phys = 0.12,
             rho0N_phys = 0.071, rho0B_phys = 1.;
// Scaling & LBM variables
const int nx = 500, ny = nx;
const int nsteps = 2000, n_out = 500;

const double Scale_length = L_phys / nx, Scale_rho = 1.,
             rho0N = rho0N_phys / Scale_rho, rho0B = rho0B_phys / Scale_rho,
             beta = beta_phys / Scale_rho;
const double tauN = 2.5, dN = (tauN - 0.5) * cs2, omegaN = 1. / tauN,
             omega1N = 1. - omegaN;
const double Scale_diffus = dN_phys / dN, dB = dB_phys / Scale_diffus,
             tauB = dB / cs2 + 0.5, omegaB = 1. / tauB, omega1B = 1. - omegaB;
const double Scale_time = pow(Scale_length, 2) / Scale_diffus;
vector<double> f1N(nx *ny *np, 0.), f2N(nx *ny *np, 0.), rhoN(nx *ny, 0.);
vector<double> f1B(nx *ny *np, 0.), f2B(nx *ny *np, 0.), rhoB(nx *ny, 0.);
vector<double> rhoD(nx *ny, 0.);
int newx, newy, id, idn, check, idk, idnk;
double bacteria, dead, nutrient, reaction, reactBact, reactNutr, bacteria0,
    coeff_reactions;
double mu_death;
const double a1 = 1. / 2400., a2 = 1. / 120.;

///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void export_ppm(std::string name, std::size_t time) {
    stringstream output_filename;
    output_filename << "out/" << name << "_" << time << ".ppm";
    ofstream of;
    of.open(output_filename.str().c_str());

    of << "P3\n";
    of << nx << " " << ny << std::endl;
    of << "255\n";

    for (size_t y = 0; y < ny; y++) {
        for (size_t x = 0; x < nx; x++) {
            int idx = y * nx + x;
            int r = (rhoN[idx] < 1e-12) ? 0 : (int)(rhoN[idx] * 8 * 255);
            int g = (rhoB[idx] < 1e-12) ? 0 : (int)(rhoB[idx] * 8 * 255);
            int b = (rhoD[idx] < 1e-12) ? 0 : (int)(rhoD[idx] * 8 * 255);

            of << r << " " << " " << g << " " << b << " ";
        }
        of << std::endl;
    }
    of.close();
}

void initial_state() {
    srand(time(NULL));
    int radius = nx / 50;
    bacteria0 = 0.;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
            id = x * ny + y;
            if (x == nx / 2 && y == ny / 2) rhoB[id] = rho0B;
            rhoN[id] = rho0N;
            bacteria0 += rhoB[id];
            for (int k = 0; k < np; k++) {
                f1N[id * np + k] = wf[k] * rhoN[id];
                f1B[id * np + k] = wf[k] * rhoB[id];
            }
        }
}

int algorithm_lattice_boltzmann(int time) {
    check = 0;
    bacteria = dead = nutrient = 0.;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
            id = x * ny + y;
            rhoN[id] = rhoB[id] = 0.;
            for (int k = 0; k < np; k++) {
                idk = id * np + k;
                rhoN[id] += f1N[idk];
                rhoB[id] += f1B[idk];
            }
            mu_death = pow((1. + rhoB[id] / a1) * (1. + rhoN[id] / a2), -1);
            coeff_reactions = 0.;
            if (rhoB[id] >= beta) coeff_reactions = rhoB[id] * rhoN[id];
			else if (rhoB[id] < beta) std::cout << "I need this" << std::endl;
            reactNutr = -coeff_reactions;
            reactBact = coeff_reactions - mu_death * rhoB[id];
            rhoD[id] += mu_death * rhoB[id];
            for (int k = 0; k < np; k++) {
                idk = id * np + k;
                f1N[idk] = omega1N * f1N[idk] + omegaN * wf[k] * rhoN[id] +
                           wf[k] * reactNutr;
                f1B[idk] = omega1B * f1B[idk] + omegaB * wf[k] * rhoB[id] +
                           wf[k] * reactBact;
                newx = x + cx[k];
                newy = y + cy[k];
                if (x == 0 || x == nx - 1) newx = (newx + nx) % nx;
                if (y == 0 || y == ny - 1) newy = (newy + ny) % ny;
                idn = newx * ny + newy;
                idnk = idn * np + k;
                f2N[idnk] = f1N[idk];
                f2B[idnk] = f1B[idk];
            }
            nutrient += rhoN[id];
            bacteria += rhoB[id];
            dead += rhoD[id];
            if (isnan(rhoN[id])) check = 1;
        }
    return check;
}

void boundaries() {
    int x;
    for (int y = 0; y < ny; y++) {
        x = 0;
        id = x * ny + y;
        rhoN[id] = 20 * rho0N;
        for (int k = 0; k < np; k++) f2N[id * np + k] = wf[k] * rhoN[id];

        x = nx - 1;
        id = x * ny + y;
        rhoN[id] = rhoN[(nx - 2) * ny + y];
        for (int k = 0; k < np; k++) f2N[id * np + k] = wf[k] * rhoN[id];
    }
}

int main(int argc, char *argv[]) {
    FILE *data = fopen("data.txt", "wt");
    system("mkdir out");
    int check_mach = 0, t;
    initial_state();
    printf("Dx = %e[m] Dt = %e[s]\n", Scale_length, Scale_time);
    printf("%e  %e %e\n", tauB, tauN, Scale_diffus);
    double scaling_factor = Scale_rho / (nx - 1) / (ny - 1);
    for (t = 0; t < nsteps; t++) {
        check_mach = algorithm_lattice_boltzmann(t);
        // boundaries();
        f1N = f2N;
        f1B = f2B;
        if (plot == true && t % n_out == 0) export_ppm("all", t);
        // if(t%n_out==0)
        printf("%d of %d. All=%e\n", t, nsteps,
               bacteria + dead);
        fprintf(data, "%e    %e        %e    %e\n", (double)(t * Scale_time),
                bacteria * scaling_factor, (bacteria + dead) * scaling_factor,
                nutrient * scaling_factor);
        if (check_mach == 1) goto labelA;
    }
labelA:
    fclose(data);
    return 0;
}