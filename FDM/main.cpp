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
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = false;
const double L = 1, dN = 1., dB = 0.05, N0 = 0.1, B0 = 1.; 
const int nx = 1000, ny = nx;
const double Dx = L/nx, Dt = 1.45e-7;
const int nsteps = 1000, n_out = 1000;
vector<double> N(nx*ny, 0.), N_old(nx*ny, 0.);
vector<double> B(nx*ny, 0.), B_old(nx*ny, 0.);
vector<double> D(nx*ny, 0.);
int id, check, xp, xm, yp, ym;
double bacteria, dead, nutrient, reaction, bacteria0, mu_death, reactBact, reactNutr, coeff_reactions, laplBact, laplNutr, h1;
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
			if(B[X*ny+Y]<1e-6)
				output_file << "0\n";
			else
				output_file << B[X*ny+Y]<< "\n";
		}

	output_file << "SCALARS Dead double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(D[X*ny+Y]<1e-6)
				output_file << "0\n";
			else
				output_file << D[X*ny+Y]<< "\n";
		}

	output_file << "SCALARS All double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	 for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(B[X*ny+Y]+D[X*ny+Y]<1e-6)
				output_file << "0\n";
			else
				output_file << B[X*ny+Y]+D[X*ny+Y]<< "\n";
		}

	output_file << "SCALARS Nutrients double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(N[X*ny+Y]<1e-6)
				output_file << "0\n";
			else
				output_file << N[X*ny+Y]<< "\n";
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
			//if( (x-nx/2)*(x-nx/2)+(y-ny/2)*(y-ny/2)<pow(radius+radius*0.001*((double)rand()/(RAND_MAX)),2) )
			if(x==nx/2 && y==ny/2)
				B[id] = B0;
			N[id] = N0;
			bacteria0 += B[id];
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
int finite_difference_algorithm(int time)
{
	check = 0;
	bacteria = dead = nutrient = 0.;
	N_old = N;
	B_old = B;
	// compute first term
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			xm = (x-1+nx)%nx;
			xp = (x+1+nx)%nx;
			ym = (y-1+ny)%ny;
			yp = (y+1+ny)%ny;
			laplNutr = (N[xm*ny+y] + N[xp*ny+y] + N[x*ny+yp] + N[x*ny+ym] - 4.*N[id])/pow(Dx,2);
			laplBact = (B[xm*ny+y] + B[xp*ny+y] + B[x*ny+yp] + B[x*ny+ym] - 4.*B[id])/pow(Dx,2);
			mu_death = pow((1.+B[id]/a1)*(1.+N[id]/a2),-1);
			h1 = Dt*dN*laplNutr - B[id]*N[id];
			N_old[id] += 0.5*h1;
			h1 = Dt*dB*laplBact + B[id]*N[id] - mu_death*B[id];
			B_old[id] += 0.5*h1;
			D[id] += 0.5*mu_death*B[id];
		}
	// compute second term
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			xm = (x-1+nx)%nx;
			xp = (x+1+nx)%nx;
			ym = (y-1+ny)%ny;
			yp = (y+1+ny)%ny;
			laplNutr = (N_old[xm*ny+y] + N_old[xp*ny+y] + N_old[x*ny+yp] + N_old[x*ny+ym] - 4.*N_old[id])/pow(Dx,2);
			laplBact = (B_old[xm*ny+y] + B_old[xp*ny+y] + B_old[x*ny+yp] + B_old[x*ny+ym] - 4.*B_old[id])/pow(Dx,2);
			mu_death = pow((1.+B_old[id]/a1)*(1.+N_old[id]/a2),-1);
			h1 = Dt*dN*laplNutr - B_old[id]*N_old[id];
			N[id] = N_old[id] + 0.5*h1;
			h1 = Dt*dB*laplBact + B_old[id]*N_old[id] - mu_death*B_old[id];
			B[id] = B_old[id] + 0.5*h1;
			D[id] += 0.5*mu_death*B_old[id];

			nutrient += N[id];
			bacteria += B[id];
			dead += D[id];
			if(isnan(N[id]))
				check = 1;
		}
    return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data = fopen("data.txt","wt");
	system("mkdir vtk_fluid");
	int check_error = 0, t;
	initial_state();
  printf("Dx = %e[m] Dt = %e[s]\n", Dx, Dt);
  for(t=0; t<nsteps; t++)
  {
    check_error = finite_difference_algorithm(t);
    if(plot_vtk==true && t%n_out==0)
			write_fluid_vtk(t);
		if(t%n_out==0)
			printf("Time %e [s] of %e [s]. Bacteria=%e\n", t*Dt, nsteps*Dt, bacteria+dead);
		fprintf(data,"%e    %e        %e    %e\n", (double)(t*Dt), bacteria/(nx-1)/(ny-1), (bacteria+dead)/(nx-1)/(ny-1), nutrient/(nx-1)/(ny-1));
		if(check_error==1)
      goto labelA;
  }
  labelA:
	fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
