///------------------------------------------------------------------------------------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
///------------------------------------------------------------------------------------------------------------------------------------------------------
const int nx = 401, ny = 401, np = 9, nsteps = 100000001, n_out = 1000000, cx[np] = {0., 1., 0., -1., 0., 1., -1., -1., 1.}, cy[np] = {0., 0., 1., 0., -1., 1., 1., -1., -1.};
const double w[np] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}, cs2 = 1./3., rho0 = 1., u_lid = 0.005, Reynolds = 7500., ni = u_lid*(ny-1)/Reynolds, tau = ni*3+.5, omega = 1./tau,  omega_vect[np] = {1., 1., 1., 1., omega, omega, 1., 1., 1.}, omega_vect1[np] = {1.-1, 1-1., 1-1., 1-1., 1-omega, 1-omega, 1-1., 1-1., 1-1.};
double f1[np][nx][ny], f2[np][nx][ny], u[nx][ny], v[nx][ny], rho[nx][ny], collision_operator[np], feq[np], k3, k4, k5, k6, k7, k8, k3_eq, k4_eq, k5_eq, k6_eq, k7_eq, k8_eq, CX, CY, CZ, U, V, R, ftemp, U2, V2, A, B, C, Pxx, Pxy, Pyy, fneq, QPI;
int newi, newj, i, j, k;
const bool plot_vtk = true;
///------------------------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	output_file << "VECTORS velocity_vector float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
			output_file << u[X][Y]/u_lid << " " << v[X][Y]/u_lid << " 0\n";

	output_file.close();
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
		{
			R = rho[i][j] = rho0;
			U = u[i][j] = 0;
            V = v[i][j] = 0.;
			C = -1.5*(U*U+V*V);
			for(k=0; k<np;k++)
			{
           		A = U*cx[k]+V*cy[k];
           		B = 4.5*A*A;
                f1[k][i][j] = f2[k][i][j] = w[k]*R*(1.+3.*A+B+C);
			}
		}
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
void my_algorithm()
{
	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++)
		{
			/// compute macroscopic variables
			R = U = V = 0.;
			for(k=0; k<np; k++)
			{
				ftemp = f1[k][i][j];
				R += ftemp;
				U += ftemp*cx[k];
				V += ftemp*cy[k];
			}
			U /= R;
			V /= R;
			u[i][j] = U;
			v[i][j] = V;
			rho[i][j] = R;
			C = -1.5*(U*U+V*V);
			/// save the post-collision state in f1 and stream in f2
			for(k=0; k<np; k++)
			{
				A = U*cx[k]+V*cy[k];
				B = 4.5*A*A;
				f1[k][i][j] = omega1*f1[k][i][j]+omega*wf[k]*R*(1.+3.*A+B+C);
				newi = i+cx[k];
				newj = j+cy[k];
				if(i==0 || i==nx-1)
					newi = (newi+nx)%nx;
				if(j==0 || j==ny-1)
					newj = (newj+ny)%ny;
				f2[k][newi][newj] = f1[k][i][j];
			}
		}
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
void boundary_conditions() /// by regularized scheme
{
	/// west and east walls
	for(j=0; j<ny; j++)
	{
    	Pxx = Pxy = Pyy = 0.;
	    V = v[0][j] = 0.;
       	U = u[0][j] = 0;
       	R = rho[0][j] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][0][j]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{		
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][0][j] = feq[k]+4.5*w[k]*QPI;
       	}
       	Pxx = Pxy = Pyy = 0.;
	    V = v[nx-1][j] = 0.;
       	U = u[nx-1][j] = 0.;
       	R = rho[nx-1][j] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][nx-1][j]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][nx-1][j] = feq[k]+4.5*w[k]*QPI;
       	}
	}
	/// south wall and north lid
	for(i=0; i<nx; i++)
	{
   		Pxx = Pxy = Pyy = 0.;
	    V = v[i][0] = 0.;
       	U = u[i][0] = 0.;
       	R = rho[i][0] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][i][0]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{		
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][i][0] = feq[k]+4.5*w[k]*QPI;
       	}
       	Pxx = Pxy = Pyy = 0.;
	    V = v[i][ny-1] = 0.;
       	U = u[i][ny-1] = u_lid;
       	R = rho[i][ny-1] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][i][ny-1]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][i][ny-1] = feq[k]+4.5*w[k]*QPI;
       	}
	 }
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
void print_velocity_profiles()
{
	FILE *ux_vertical = fopen("ux_vertical_mid_section.txt","wt"); 
	for(int j=0; j<ny; j++)
		fprintf(ux_vertical, "%d  %lf\n", j, u[nx/2][j]/u_lid);
	fclose(ux_vertical);
	
	FILE *uy_horizontal = fopen("uy_horizontal_mid_section.txt","wt"); 
	for(int i=0; i<nx; i++)
		fprintf(uy_horizontal, "%d  %lf\n", i, v[i][ny/2]/u_lid);
	fclose(uy_horizontal);
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	system("mkdir vtk_fluid");
	initial_state();
	for(int t=0; t<nsteps; t++)
    {
		my_algorithm(); // compute macroscopic variables, collide and stream
		boundary_conditions(); // apply boundary conditions 
		memcpy(f1, f2, sizeof(f2)); // swap f2 in f1
		if(plot_vtk==true && t%n_out==0)
		{
			printf("Iteration %d of %d.\n", t, nsteps);
			write_fluid_vtk(t);
		}
    }
	print_velocity_profiles();
    return 0;
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------------------------
