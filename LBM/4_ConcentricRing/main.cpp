///----------------------------------------------------------------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
// https://www.sciencedirect.com/science/article/pii/S0378437198003458
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = true;
const double cs2 = 1./3.;
const int np = 5;
const vector<int> cx = {0, 1, 0, -1, 0},
									cy = {0, 0, 1, 0, -1},
								 opp = {0, 3, 4, 1, 2};
const vector<double> wf = {1/3., 1/6., 1/6., 1/6., 1/6.};
/*
const int np = 9;
const vector<int> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1},
								  cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
*/
const double L_phys = 1, beta_phys = 0., dN_phys = 1., dB_phys = 0.05, rho0N_phys = 0.1, rho0B_phys = 1.; 
// Scaling & LBM variables
const int nx = 2000, ny = nx;
const double Scale_length = L_phys/nx, Scale_rho = 1., rho0N = rho0N_phys/Scale_rho, rho0B = rho0B_phys/Scale_rho, beta = beta_phys/Scale_rho;
const double tauN = 2.5, dN = (tauN-0.5)*cs2, omegaN = 1./tauN, omega1N = 1.-omegaN;
const double Scale_diffus = dN_phys/dN, dB = dB_phys/Scale_diffus, tauB = dB/cs2+0.5, omegaB = 1./tauB, omega1B = 1.-omegaB;
const double Scale_time = pow(Scale_length,2)/Scale_diffus;
const int nsteps = 100000, n_out = 1000;
vector<double> f1N(nx*ny*np, 0.), f2N(nx*ny*np, 0.), rhoN(nx*ny, 0.);
vector<double> f1B(nx*ny*np, 0.), f2B(nx*ny*np, 0.), rhoB(nx*ny, 0.);
vector<double> rhoD(nx*ny, 0.);
int newx, newy, id, idn, check, idk, idnk;
double bacteria, dead, nutrient, reaction, reactBact, reactNutr, bacteria0, coeff_reactions;
double mu_death;
const double a1 = 1./2400., a2 = 1./120.;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " double\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " double\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " double\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	output_file << "SCALARS Bacteria double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoB[X*ny+Y]<1e-12)
				rhoB[X*ny+Y] = 0.;
			output_file << rhoB[X*ny+Y]*Scale_rho << "\n";
		}

	output_file << "SCALARS Dead double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoD[X*ny+Y]<1e-12)
				rhoD[X*ny+Y] = 0.;
			output_file << rhoD[X*ny+Y]*Scale_rho << "\n";
		}

	output_file << "SCALARS All double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoB[X*ny+Y]+rhoD[X*ny+Y]<1e-12)
				rhoB[X*ny+Y] = rhoD[X*ny+Y] = 0.;
			output_file << (rhoB[X*ny+Y]+rhoD[X*ny+Y])*Scale_rho << "\n";
		}

	output_file << "SCALARS Nutrients double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoN[X*ny+Y]<1e-12)
				rhoN[X*ny+Y] = 0.;
			output_file << rhoN[X*ny+Y]*Scale_rho << "\n";
		}
	
	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	srand(time(NULL));
	int radius = nx/50;
	bacteria0 = 0.;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			if( x==nx/2 && y==ny/2)
				rhoB[id] = rho0B;
			rhoN[id] = rho0N;
			bacteria0 += rhoB[id];
			for(int k=0; k<np; k++)
			{
				f1N[id*np+k] = wf[k]*rhoN[id];
				f1B[id*np+k] = wf[k]*rhoB[id];
			}
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algorithm_lattice_boltzmann(int time)
{
	check = 0;
	bacteria = dead = nutrient = 0.;
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			rhoN[id] = rhoB[id] = 0.;
			for(int k=0; k<np; k++)
			{
				idk = id*np+k;
				rhoN[id] += f1N[idk];
				rhoB[id] += f1B[idk];
			}
			mu_death = pow((1.+rhoB[id]/a1)*(1.+rhoN[id]/a2),-1);
			coeff_reactions = 0.;
			if(rhoB[id]>=beta)
				coeff_reactions = rhoB[id]*rhoN[id];
			reactNutr = -coeff_reactions;
			reactBact = coeff_reactions - mu_death*rhoB[id];
			rhoD[id] += mu_death*rhoB[id];
			for(int k=0; k<np; k++)
			{
				idk = id*np+k;
				f1N[idk] = omega1N*f1N[idk] + omegaN*wf[k]*rhoN[id] + wf[k]*reactNutr;
				f1B[idk] = omega1B*f1B[idk] + omegaB*wf[k]*rhoB[id] + wf[k]*reactBact;
				newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
				idnk = idn*np+k;
				f2N[idnk] = f1N[idk]; 
				f2B[idnk] = f1B[idk];
			}
			nutrient += rhoN[id];
			bacteria += rhoB[id];
			dead += rhoD[id];
			if(isnan(rhoN[id]))
				check = 1;
		}
    return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
void boundaries()
{
	int x;
	for(int y=0; y<ny; y++)
	{
		x = 0;
		id = x*ny+y;
		rhoN[id] = 20*rho0N;
		for(int k=0; k<np; k++)
			f2N[id*np+k] = wf[k]*rhoN[id];
		
		x = nx-1;
		id = x*ny+y;
		rhoN[id] = rhoN[(nx-2)*ny+y];
		for(int k=0; k<np; k++)
			f2N[id*np+k] = wf[k]*rhoN[id];
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data = fopen("data.txt","wt");
	system("mkdir vtk_fluid");
	int check_mach = 0, t;
	initial_state();
	printf("Dx = %e[m] Dt = %e[s]\n", Scale_length, Scale_time);
  printf("%e  %e %e\n", tauB, tauN, Scale_diffus);
  double scaling_factor = Scale_rho/(nx-1)/(ny-1);
  for(t=0; t<nsteps; t++)
  {
    check_mach = algorithm_lattice_boltzmann(t);
    //boundaries();
    f1N = f2N;
    f1B = f2B;
		if(plot_vtk==true && t%n_out==0)
			write_fluid_vtk(t);
		//if(t%n_out==0)
			printf("%e of %e. All=%e\n", t*Scale_time, nsteps*Scale_time, bacteria+dead);
		fprintf(data,"%e    %e        %e    %e\n", (double)(t*Scale_time), bacteria*scaling_factor, (bacteria+dead)*scaling_factor, nutrient*scaling_factor);
		if(check_mach==1)
      goto labelA;
  }
  labelA:
	fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
