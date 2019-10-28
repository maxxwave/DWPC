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

#include "../hdr/storage.hdr"
#include "../hdr/calculate.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/rc.hdr"

namespace reservoir{
	// Define the number of neurons alias nodes
	int no_nodes=24;
	//define the time for a 
	// In this variable we set the time need for a single discrete input
	double tau=4.3e-7; //s
	// In this variable we se how long is the time applied for a single node
	double theta=tau/no_nodes;

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
	double V_min=10;
	double V_max=40;
	double oscillation_response(){
       		// In this loop we apply a sequence of input fields from 	
		for (stor::V0=V_min; stor::V0<=V_max; stor::V0 +=V_max*0.25){
			//we loop over the nodes
			for (int i=1; i<=no_nodes;i++){
				// In this stage we apply the mask to the input signal
				double mask_rand=mask();
				stor::V = stor::V0*mask_rand;
				for (i=0; i<=theta; i+=integrate::Dt)
					integrate::runge_kutta(time);
					std::cout<<stor::V<<"\t"<<stor::V0<<"\t"<<
                        		stor::x_dw<<"\t"<<time*1e9<<std::endl;
			}

		}
	}


}
