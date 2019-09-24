// Author: Razvan Ababei
// date: 15.09.2019
// University of Sheffield

#include <iostream>
#include <cmath>
#include <fstream> 

#include "../hdr/calculate.hdr"
#include "../hdr/storage.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/program.hdr"

namespace programs{
	// in this function we plot the bifurcation diagram of domain wall propagation 
	double bifurcation(){
		std::ofstream output;
		output.open("bifurcation.data");
		// In this variable we store the no of integration steps performed within one period of the applied field
		double no = 1/(stor::omega*integrate::Dt);
		std::cout<<no<<std::endl;
		
		// the applied field amplitude is varied from 1000A/M to 5000A/M with a step of 4A/m
	       	for(int i=0;i<1000;i++){
			stor::V0=i*4+1000;

			// This is a local variable dedicated to control the time simulation
			double time=0.0;
			//perform some equilibration steps
			for (int k=0; k<=500000; k++){
				time += integrate::Dt;
				calculate::Zeeman(time);
				integrate::euler();
			}//end of equilibration

			// In this loop we perform a numebr of integration steps equal to a period
			for (int l=0; l<100; l++){
				for (int t=0; t<=2000000;t++){
					time += integrate::Dt;
				       	calculate::Zeeman(time);
					integrate::euler();		
				}
				// print the value of the DW position after each loop (period)
				output << "H="<<"\t"<<stor::V<<"\t"<<"X="<<"\t"<<stor::x_dw<<"\t"<<"V="<<"\t"<<stor::vx<<std::endl;

			}// end of loop

		
		}	
		
		output.close();
	}//end of bifurcation function
	
}// end of namespace
