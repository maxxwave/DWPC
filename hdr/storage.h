#ifndef __STORAGE_H__
#define __STORAGE_H__

#include <iostream>
#include <vector>
#include <string>

#include "../hdr/arrays.h"


const double Pi =  3.1415926535897932384626;
const double Pi_sqr =  Pi*Pi;

namespace stor{

    extern std::string program;
	extern double L;
	extern double Lz, Ly;
	extern double Ms, mu0;
	extern double muMs;
	extern double Nx, Ny, Nz;
	extern double alpha;
	extern double A;
	extern double xl,xr;
	extern double V, V0, freq, omega;
    extern double H_const;
	extern double cell_size;
	extern int Nwires;
	extern std::vector<double> V0_mdw;
	extern std::vector<double> x_coord;
	extern std::vector<double> phi_coord;
	extern double T_sim;
	extern double gamma;
	//extern std::vector<double> H_demag;
	extern double H_demag;
	extern double Delta;
	extern double x_dw;
	extern double phi_dw;
	extern double Ex, dEx;
	extern double vx, phi_dt;
	extern double Dw_size;
	extern double u, mu_B, e_el, beta, P, j_dens;
	extern double rij;
	extern std::vector<double> H_DW;
	extern std::vector<double> H_dd;
	extern std::vector<double> j_dens_dw;
	extern std::vector<double> u_dw;
	extern double my;
	extern void initialize();
	extern double A0, A1, A2, A3, A4, A5, A6, A7, A8;
    extern bool use_DW_coupling;

    extern bool use_edge;
    extern double H_edge_max;
    extern double edge_scale;
    extern array_t<2,double> H_edge;
    extern array_t<2,double> E_edge;
} //end of namespace
#endif
