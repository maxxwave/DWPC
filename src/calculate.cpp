// Author: Razvan Ababei
// University of Sheffield
// date:20th of August, 2019
//
//
// This file is dedicated to create the chain including the potential for strip
// This model includes two anti-notches where the potential follows a form given by:
//
//


#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <string>

#include "../hdr/storage.h"
#include "../hdr/calculate.h"
#include "../hdr/arrays.h"
#include "../hdr/euler_integrator.h"

namespace calculate{
	// Defining some prefactors where we incorporate the constants in order to not be called each time in the loop
	double prefac1 = (-stor::alpha*stor::gamma)/((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
	double prefac2 = stor::mu0*stor::gamma*stor::H_demag/2.0; //(2.0+2.0*stor::alpha*stor::alpha);
	double prefac3 = -stor::gamma/((1+ stor::alpha*stor::alpha)*2*stor::Ms*stor::Lz*stor::Ly);
	double prefac4 = -(stor::gamma*stor::alpha*stor::mu0*stor::H_demag)/(2+2*stor::alpha*stor::alpha);
	double zeeman_prefac1 = stor::gamma*stor::mu0*stor::alpha/(stor::alpha*stor::alpha+1.0);
	double zeeman_prefac2 = stor::gamma*stor::mu0/(1.0 + stor::alpha*stor::alpha);
	const double one_rad=Pi/180.0;
	double current_prefac=(1+2*stor::alpha*stor::beta-stor::alpha*stor::alpha);
	// Calculating the number of cells for a given L and cell_size
	int a = int(stor::L/stor::cell_size);
	int N=2*a+1;
	double kb=1.38064e-23;
	// some parameters from origin fit
	// The expression is F(X)= A0 + A1*X**2 + A2*X**2 + A3*X**3 + A4*X**4 + ... + A8*X**8
	// These values correspond to Py antinotches
	/*double a0 = 2.395e-20;
	double a1 = 6.741e-15;
	double a2 = -1.29e-7;
	double a3 = -2.94244;
	double a4 = -2.567e8;
	double a5 = 2.91623e14;
	double a6 = 6.77e21;
	double a7 = -7.728e27;
	double a8 = 2.6771e35;*/
	//These values correspond to Ni antinotches
	double a0 =0;// 2.21117e-21;
	double a1 =0;// -3.8271e-15;
	double a2 =0;// -1.2866e-6;
	double a3 =0;//  0.61164;
	double a4 =0;//  1.632e8;
	double a5 = 0.0;
	double a6 = 0.0;
	double a7 = 0.0;
	double a8 = 0.0;

    // function which calculate the potential energy depending on x
    double update_energy_antinotches(double x){

        // this is the potential for notches
        // need to find out who is a and b ???
        //stor::Ex = stor::A0 + stor::A1*x + stor::A2*x*x + stor::A3*pow(x,3) + stor::A4*pow(x,4) + stor::A5*pow(x,5) + stor::A6*pow(x,6) + stor::A7*pow(x,7)+ stor::A8*pow(x,8);

        // this is the analytical derivative of the potential dE/dx
        //stor::dEx=stor::A1+ 2*stor::A2*x + 3*stor::A3*x*x + 4*stor::A4*pow(x,3) + 5*stor::A5*pow(x,4) + 6*stor::A6*pow(x,5) + 7*stor::A7*pow(x,6) + 8*stor::A8*pow(x,7);
        stor::dEx= 2*x*(stor::A2 + x*( 2*stor::A4*x + 3*stor::A6*x*x*x) );

        return stor::dEx;
    }// end of function


    double compute_Hpin(double &x)
    {
        return 2*x*(stor::A2 + x*( 2*stor::A4*x + 3*stor::A6*x*x*x) );
    }

    double Vp_2deriv( const double x) {
        return 2*stor::A2 + 12*stor::A4*x*x;
    }

    // in this function we will calculate the pinning energy for anti-notches
    /*double update_energy_notches(double x){
      stor::Ex=U0-U1*(exp(-(double x + stor::xl)*(double x + stor::xl)/(stor::L*stor::L) )
      -(exp(-(double x + stor::xl)*(double x + stor::xr)/(stor::L*stor::L) )
      stor::dEx=-2*U1/(stor::L*stor::L)*(x+stor::xl)*exp(-(x+stor::xl)*(x+stor::xl)/(stor::L*stor::L)) +
      +2*U1/(stor::L*stor::L)*(x+stor::xr)*exp(-(x+stor::xr)*(x+stor::xr)/(stor::L*stor::L))

      }// end of function
      */
    double update_energy();

    double K_eff( double phi) {
        return stor::muMs*stor::Ms* sin(phi)*sin(phi) + stor::muMs*stor::H_demag;
    }

    double dK_eff(double phi) {
        return 2.0*stor::muMs * stor::Ms *sin(phi)*cos(phi);
    }

    double DW(double phi) {
        return Pi * sqrt( 2.0 * stor::A / K_eff(phi));
    }
    //In this function we implement the stochastic term
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);
    unsigned random_seed = 1234;
    void seed_rng(const unsigned seed) {
        generator.seed(seed);
    }

    double Normal() {
        return distribution(generator);
    }

    void Normal(std::vector<double> &x) {
        for ( int i=0; i < x.size(); i++)
            x[i] = distribution(generator);
    }

    double noise(double T, double DW){
        double rand=distribution(generator);
	//std::cout<<"temperatura setata=  "<<rand<<std::endl;
    	double noise= sqrt((2*stor::alpha*kb*T*stor::gamma/(stor::Ms*stor::Lz*stor::Ly*DW*integrate::Dt))) * rand; //1/integrate::dt
	return noise;
    }

    // function which calculate the Domain wall width
    double calculate_DW(double phi){
        stor::Dw_size = Pi * sqrt(2*stor::A / K_eff(phi));
        // Pivano form of DW
        //stor::Dw_size=sqrt(2*stor::A/(stor::mu0*stor::Ms*stor::Ms*(stor::Ny*sin(phi)*sin(phi) + stor::Nz*cos(phi)*cos(phi)))); // Matt form
        return stor::Dw_size;
    }// end of function calculate_DW

    double DW_gradient( double phi, double Delta) {
        return - Pi_sqr * (stor::A / ( DW(phi) * K_eff(phi)*K_eff(phi))) * dK_eff(phi);
    }

    // In this routine we calculate the Zeeman field taking into account the frequency of the field
    double Zeeman(double time){

        stor::V = stor::V0*sin(stor::omega*time);

		//this equation can be used for benchmark1 program
		//stor::V=stor::V0*cos(stor::omega*time);
		return stor::V;
	}

    double Zeeman(double time, const int i){
        //std::cout<<stor::V0_mdw[0]<<"\t"<<stor::V0_mdw[1]<<std::endl;
        return stor::V0_mdw[i]*sin(stor::omega*time);
    }

    // we define two function for speed and angular speed
    double phi_t(double dEx, double phi_rk, double H){
        return prefac3*dEx+prefac4*sin(2*phi_rk)+ zeeman_prefac2*H;
    }
    double x_t(double DWs, double phi, double phi_t){
        return prefac2*sin(2*phi)*DWs + stor::alpha*DWs*phi_t;
    }

    double current( double time, double j_dens ){
		//std::cout<<stor::P<<"\t"<<stor::mu_B<<"\t"<<stor::j_dens<<std::endl; //need to add a time function
		//return stor::j_dens*stor::P*stor::mu_B/(stor::e_el*stor::Ms)*exp(-((time*time+25e-18-2*time*5e-9)/1e-18)); //need to add a time function
		//if(time>=5e-9){
		return stor::j_dens*j_dens*stor::P*stor::mu_B/(stor::e_el*stor::Ms);
	}

	// in this routine we calculate the DWs coupling
	double DW_coupling( std::vector <double> &X_DW ){
		if(stor::Nwires>1){
		// the monopole interaction
		stor::H_DW.resize(stor::Nwires);
		// the dipolar interaction
		stor::H_dd.resize(stor::Nwires);

		double S=stor::Ly*stor::Lz;
		double rijd=(stor::rij+stor::Ly)*(stor::rij+stor::Ly);
		double prefac_dd=0.375*stor::my*stor::my/(Pi*stor::Ms*S);

		// loop over the wires
		for (int i=0; i<X_DW.size(); i++){
			double r=sqrt((X_DW[i]-X_DW[i+1])*(X_DW[i]-X_DW[i+1]) + rijd);
			double r3=r*r*r;
			double r7=r*r*r*r*r*r*r;

			//boundary wires
			double r_sec=sqrt((X_DW[0]-X_DW[1])*(X_DW[0]-X_DW[1]) + rijd);
			double r_sec3=r_sec*r_sec*r_sec;
			double r_sec7=r_sec*r_sec*r_sec*r_sec*r_sec*r_sec*r_sec;

			stor::H_dd[0]= prefac_dd*(X_DW[0]-X_DW[1])*((X_DW[0]-X_DW[1]) *(X_DW[0]-X_DW[1])-4*rijd)/r_sec7;

			stor::H_DW[0] = -stor::Ms*S*(X_DW[0]-X_DW[1])/(2*Pi*r_sec3);

			double r_prim=sqrt((X_DW[X_DW.size()-1]-X_DW[X_DW.size()-2])*(X_DW[X_DW.size()-1]-X_DW[X_DW.size()-2]) + rijd);
			double r_prim3=r_prim*r_prim*r_prim;
			double r_prim7=r_prim*r_prim*r_prim*r_prim*r_prim*r_prim*r_prim;

			stor::H_dd[X_DW.size()-1]= prefac_dd*(X_DW[X_DW.size()-1]-X_DW[X_DW.size()-2])*
				((X_DW[X_DW.size()-1]-X_DW[X_DW.size()-2]) * (X_DW[X_DW.size()-1]-X_DW[X_DW.size()-2])-4*rijd)/r_prim7;

			stor::H_DW[X_DW.size()-1] = -stor::Ms*S*(X_DW[X_DW.size()-1]-X_DW[X_DW.size()-2])/(2*Pi*r_prim3);

			if((i!=0)&&(i!=(X_DW.size()-1))){
			// We assume the NN interaction only
			stor::H_DW[i] = -stor::Ms*S*(X_DW[i]-X_DW[i+1])/(2*Pi*r3)
					-stor::Ms*S*(X_DW[i]-X_DW[i-1])/(2*Pi*r3);

			stor::H_dd[i]= prefac_dd*(X_DW[i]-X_DW[i+1])*((X_DW[i]-X_DW[i+1])*(X_DW[i]-X_DW[i+1])-4*rijd)/r7
				     + prefac_dd*(X_DW[i]-X_DW[i-1])*((X_DW[i]-X_DW[i-1])*(X_DW[i]-X_DW[i-1])-4*rijd)/r7;
			}

		}
	}
	return 0;
	}

    void gradient ( double &dx, double &dphi, double x, double phi, const double time)
    {
        //const double dEx = update_energy_antinotches(x);
        //const double dEx = compute_Hpin(x);
        double dEx = update_energy_antinotches(x);
        double H = Zeeman(time);
        double DWs = calculate_DW(phi);
        double n_x = noise(stor::T_sim, DWs);
        double n_phi= noise(stor::T_sim, DWs);
        double u=current(time,stor::j_dens);
        dphi = prefac3*dEx + prefac4*sin(2*phi)
            + zeeman_prefac2*H //+ (n_phi + stor::alpha*n_x)/(1+stor::alpha*stor::alpha)
            + (stor::beta-stor::alpha)*u/DWs;
        dx = prefac2*sin(2*phi)*DWs + stor::alpha*DWs*dphi
            + u;
        //std::cout<<u<<"\t"<<time<<std::endl;
        //std::cout<<"nx=   "<<n_x<<"\t"<<n_phi<<"\t"<<prefac3*dEx<<std::endl;
        //double d = 0.2, g = 0.3, a = 1, b = -1, w=1;
        //dx = phi;
        //dphi = -d*phi - b*x -a*x*x*x + g*cos(w*time);

    }



    void noise_gradient(double &gx, double &gp, double x, double phi, double n_x, double n_phi)
    {
        double DWs = calculate_DW(phi);
        double gamma_p = stor::gamma/(1+stor::alpha*stor::alpha);
        double sigma = sqrt((stor::alpha*calculate::kb*stor::T_sim*gamma_p/(stor::Ms*stor::Lz*stor::Ly*DWs*integrate::Dt)));
        //gp = sigma*(n_phi - stor::alpha*n_x)/(1+stor::alpha*stor::alpha);
        //gx = sigma*(n_x + stor::alpha*n_phi)*DWs/(1+stor::alpha*stor::alpha);
        gp = sigma*n_phi;
        gx = sigma*n_x*DWs;
    }



    void gradient ( std::vector<double> &dx, std::vector<double> &dphi, std::vector<double> &x, std::vector<double> &phi, const double time)
    {
        // calculate the DW coupling
        //DW_coupling(integrate::multi_dw::x_p);
        if(stor::use_DW_coupling) DW_coupling(x);
        //std::cout<<stor::x_coord[0]<<std::endl;
	//DW_coupling(x);

        for ( int i = 0; i < x.size(); i++) {
            double dEx = update_energy_antinotches(x[i]);
            double H = Zeeman(time, i) - stor::H_DW[i];
            double DWs = calculate_DW(phi[i]);
<<<<<<< HEAD
            
	    double n_x = 0.0; // noise(stor::T_sim, DWs);
            double n_phi=0.0; // noise(stor::T_sim, DWs);
=======
            double n_x =0.0;  // noise(stor::T_sim, DWs);
            double n_phi=0.0; //  noise(stor::T_sim, DWs);
>>>>>>> bfcd5d580d58c958e44044a4b890c02297d5e505

            stor::u_dw[i]=current(time, stor::j_dens_dw[i]);

            dphi[i] = prefac3*dEx + prefac4*sin(2*phi[i]) + zeeman_prefac2*H
                + (stor::beta-stor::alpha)
		*stor::u_dw[i]/DWs;
            dx[i] = prefac2*sin(2*phi[i])*DWs + stor::alpha*DWs*dphi[i]
                + stor::u_dw[i];
        }
    }

    void noise_gradient(std::vector<double> &gx, std::vector<double> &gp, std::vector<double> &x, std::vector<double> &phi, std::vector<double> &n_x, std::vector<double> &n_phi)
    {
        for ( int i = 0; i < gx.size(); i++){
            double DWs = calculate_DW(phi[i]);
            double gamma_p = stor::gamma/(1+stor::alpha*stor::alpha);
            double sigma = sqrt((stor::alpha*calculate::kb*stor::T_sim*gamma_p/(stor::Ms*stor::Lz*stor::Ly*DWs*integrate::Dt)));
            //gp = sigma*(n_phi - stor::alpha*n_x)/(1+stor::alpha*stor::alpha);
            //gx = sigma*(n_x + stor::alpha*n_phi)*DWs/(1+stor::alpha*stor::alpha);
            gp[i] = sigma*n_phi[i];
            gx[i] = sigma*n_x[i]*DWs;
        }
    }


    // Routine to compute the gradient of the differential equation (Jacobian)
    // J =  ( dfx/dx   ,   dfx/dphi )
    //      ( dfphi/dx ,   dfphi/dphi)
    void Jacobian( array_t<2,double> &J, double x, double phi, const double time)
    {
        double dEx = update_energy_antinotches(x);
        double H = Zeeman(time);
        double DWs = calculate_DW(phi);
        double dphi = prefac3*dEx + prefac4*sin(2*phi) + zeeman_prefac2*H;

        double dDelta = DW_gradient(phi, DWs);

        J(1,1) = 2.0*prefac4*cos(2.0*phi);  // dfphi/dphi
        J(0,0) = stor::alpha*DWs*prefac3*Vp_2deriv(x); //dfx/dx

        J(0,1) = prefac2*(2*DWs*cos(2*phi) + sin(2*phi)*dDelta) + stor::alpha*dDelta*dphi + stor::alpha*DWs*J(1,1);
        J(1,0) = prefac3*Vp_2deriv(x);


        //double d = 0.2, g = 0.3, a = 1, b = -1, w=1;
        //J(0,0) = 0.0;
        //J(0,1) = 1.0;
        //J(1,0) = -b - 3*a*x*x;
        //J(1,1) = -d;

    }



}//end of namespace


