// This file is dedicated to initialize the program variables from a file allowing the user to set preffered parameters
// This file takes part from DWPC code
// 
// 
// (c) author: Razvan Ababei
// University of Sheffield
// date: 11.08.2019
//
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include "../hdr/storage.hdr"
#include "../hdr/euler_integrator.hdr"

#define Pi  3.1415926535897932384626433832795028841971693993751058209749445923078164062
namespace stor{
	double initialize(){
		std::ifstream input ("input");
  		std::string line;
		std::string test1, test2, test3, test4, test5, test6, test7, test8, test9;
		std::string::size_type sz;
		if (input.is_open()){
			while ( getline (input,line) ){
				test1="Ms=";
     				if( (strstr(line.c_str(),test1.c_str()))){
					stor::Ms=std::stod(line.substr(3));
				}
				test2="L=";
     				if( (strstr(line.c_str(),test2.c_str()))){
					stor::L=std::stod(line.substr(2));
				}
				test3="Ly=";
     				if((strstr(line.c_str(),test3.c_str()))){
					stor::Ly=std::stod(line.substr(3));
				}
				test4="Lz=";
     				if((strstr(line.c_str(),test4.c_str()))){
					stor::Lz=std::stod(line.substr(3));
				}
				test5="Aex=";
     				if((strstr(line.c_str(),test5.c_str()))){
					stor::A=std::stod(line.substr(4));
				}
				test6="alpha=";
     				if((strstr(line.c_str(),test6.c_str()))){
					stor::alpha=std::stod(line.substr(6));
				}
				test7="H=";
     				if((strstr(line.c_str(),test7.c_str()))){
					stor::V0=std::stod(line.substr(2));
				}
				test8="f=";
     				if((strstr(line.c_str(),test8.c_str()))){
					stor::freq=std::stod(line.substr(2));
					stor::omega=2*Pi*stor::freq;
				}
				test9="Dt=";
     				if((strstr(line.c_str(),test9.c_str()))){
					integrate::Dt=std::stod(line.substr(3));
				}
				if( strstr(line.c_str(), "totaltime=")){
					integrate::totaltime = std::stod(line.substr(10));
				}
				if( strstr(line.c_str(), "out_time=")){
					integrate::out_time = std::stod(line.substr(9));
				}
				if( strstr(line.c_str(), "integrator=")){
					integrate::scheme = line.substr(11);
				}

				// and few more others
     					
    			}//end of while loop
			
			std::cout<<"The program has been initialized with following parameters:"<<std::endl;
			std::cout<<"Saturation, Ms = "<<stor::Ms<<"A/m"<<std::endl;			
			std::cout<<"Length of the strip = "<<stor::L<<"m"<<std::endl;			
			std::cout<<"Width = "<<stor::Ly<<"m"<<std::endl;			
			std::cout<<"Thickness = "<<stor::Lz<<"m"<<std::endl;			
			std::cout<<"Exchange stiffness, A = "<<stor::A<<"J/m"<<std::endl;			
			std::cout<<"Gilbert damping, alpha = "<<stor::alpha<<""<<std::endl;			
			std::cout<<"Amplitude of oscillating field, H0 = "<<stor::V0<<"A/m"<<std::endl;			
			std::cout<<"Frequency of the field, omega = " <<stor::omega<<"Hz"<<std::endl;			   
			std::cout<<"Integration time step, Dt = "<<integrate::Dt <<"s"<<std::endl;		
			std::cout<<"Initialization completed!"<<std::endl;
			std::cout<<"=====================================================================<"<<std::endl;	
		
		}
		else{
			std::cout<<"The program couldn't open the input file! Check the spelling of the name. This must be named 'input'"<<std::endl;
		}
		input.close();
	
	}
}//end of namespace
