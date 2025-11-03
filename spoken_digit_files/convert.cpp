#include <iostream> 
#include <fstream> 
#include "../hdr/arrays.h"

int main(){

	std::ifstream file;
	file.open("../spoken_digit_files/Theo.dat");
	array_t<2,double> sp_sig;
	sp_sig.assign(500,1025,0);
	for (int i=0; i<500; i++){
		for (int j =0; j<1025; j++){
			file>>sp_sig(i,j);
			std::cout<<sp_sig(i,j)<<std::endl;

		}
	}	
	file.close();
	return 0;
}
