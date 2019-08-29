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
	/*for (int i=0; i<stor::E_x.size(); i++){
		std::cout<< stor::E_x[i]<<"\t"<<stor::x_coord[i]<<std::endl;
	}*/

	// perform some integrations
	for (int i=0; i<100000; i++){
		for(int j=0; j<10000; j++){
			integrate::euler();
		}
		
		std::cout<<stor::x_dw<<"\t"<<stor::phi_dw<<"\t"<<i*10000*integrate::Dt<<std::endl;
	}	
return 1;	

}// end of program 
