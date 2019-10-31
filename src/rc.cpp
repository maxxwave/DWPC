// This file is dedicated to simulate the reservoir computing using a nano-strip DW oscillator 
// Author: Razvan Ababei
// University of Sheffield
// date: 28th Oct, 2019
// This file is a part of DWPC code dedicated to simulate the domain-wall propagation for a range of field inputs aiming to reproduce a neuronal network
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include "../hdr/storage.hdr"
#include "../hdr/calculate.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/rc.hdr"

namespace reservoir{

	// Define a vector to store the amplitudes of the field
	// 
	std::vector <double> H0{5,15,20,15,5,-10,-20,-10,20,20,20,20,-20,-20,-20,-20}; //in Oe

	// Define the number of neurons alias nodes
	int no_nodes=24;
	//define the time for a 
	// In this variable we set the time need for a single discrete input
	double tau=4.3e-6; //s
	// In this variable we se how long is the time applied for a single node
	double theta=tau/no_nodes;

	// we define a variable to store the no of steps needed to be performed on each node
	int long no_steps_per_node=0;

	// In this function we create a mask for the input signal 
	// This mask has a number of neurons defined above with specific weights w1, w2, w3 ...
	// This soubroutine will generate a random random either +1 or -1 only
	double mask(){
		// local random no
		double n=0;
		n=rand()%2;
		if(n==1){
			n=1;}
		else
		{n=-1;}

		return n;
	}	
	// In this routine we get the oscillator response x_i(t), where i is the sequential node
	// i=0..24
	double time=0.0;
	double oscillation_response(){
		// define a output file to store the data
		std::ofstream outputfile;
        	outputfile.open("reservoir.data");


		// we calculate the no of steps needed to be performed per node
		no_steps_per_node=std::round(theta / integrate::Dt);

       		// In this loop we apply a sequence of input fields from 	
		for (int t=0; t<H0.size(); t++){
			stor::V0=H0[t]*80;

			//we loop over the nodes
			for (int i=1; i<=no_nodes;i++){

				// In this stage we apply the mask to the input signal
				double mask_rand=mask();
				stor::V0 *= mask_rand;
				//std::cout<<no_steps_per_node<<"\t"<<theta<<"\t"<<integrate::Dt<<std::endl;
				//std::cout<<stor::V0<<std::endl;
				
				// In this loop we average over a time=theta
				for (int j=0; j<=no_steps_per_node; j++){
					integrate::runge_kutta(time);
				//	std::cout<<stor::V<<"\t"<<stor::V0<<"\t"<<
                        	//	stor::x_dw<<"\t"<<time*1e9<<std::endl;
				}
				
				outputfile << std::setprecision(4)
					   << stor::V<< "\t"
					   << stor::V0<<"\t"
					   << stor::x_dw<<"\t"
					   << time*1e9<<std::endl;
			}
		}

	outputfile.close();
	}
}
