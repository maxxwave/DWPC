// Author: Razvan Ababei
// University of Sheffield
// date 20th of August, 2019
//
//
// This file is dedicated to create the main part of the Domain-Wall Propagation Code (DWPC)
//
//
//

#include <iostream>
#include <cstdlib>

// include the headers
#include "../hdr/storage.hdr"
#include "../hdr/create.hdr"
#include "../hdr/euler_integrator.hdr"
//#include "euler_integrator.cpp"
//int main(int argc, char* argv[])

int main(){
	// execute the creation function
	create::create();
	// execute the initialization function
	create::initialize();

	// perform some integrations
	for (long int i=0; i<100000; i++){
		for(long int j=0; j<10000; j++){
			create::Zeeman(i*j*integrate::Dt);
			integrate::euler();
		}
		
	std::cout<<"x_dw="<<"\t"<<stor::x_dw<<"\t"
		<<"phi="<<"\t"<<stor::phi_dw<<"\t"
		<<"vx="<<"\t"<<stor::vx<<"\t"
		<<"t="<<"\t"<<10000*i*integrate::Dt<<std::endl;
	}	
return 1;	

}// end of program 
