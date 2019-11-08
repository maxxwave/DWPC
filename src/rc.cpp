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
#include <random>

#include "../hdr/storage.hdr"
#include "../hdr/calculate.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/rc.hdr"
#include "../hdr/input_map.h"

namespace reservoir{

    std::default_random_engine rng;
    std::uniform_real_distribution<double> rng_dist(0.0,1.0);
    std::uniform_int_distribution<int> rng_int_dist(0,1);

	std::ofstream outputfile;

	// Define a vector to store the amplitudes of the field
    //Hc defines the centre field amplitude
    double Hc = 1000;
    // dH defines the width of the field amplitude
    double dH = 100;
	std::vector <double> H0{0.25,0.75,1,0.75,0.25,-0.5,-1,-0.5,1,1,1,1,-1,-1,-1,-1}; //in Oe

    //define a target vector for the input
	std::vector <double> t_p{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
	// Define the number of neurons alias nodes
	int no_nodes=24;
	//define the time for a
	// In this variable we set the time need for a single discrete input
	double tau=12e-8; //s
	// In this variable we se how long is the time applied for a single node
	double theta=tau/no_nodes;
	// In this vector we store the arrays of outputs
	std::vector <double> s_x;
	// define an array for mask values
	std::vector <double> mask_array(no_nodes,0);

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

    // in this function we generate random number for the mask which can be either -1 or +1
	double mask_values(){
        mask_array.assign( no_nodes, 0);
		for (int j=0; j<no_nodes; j++){
			mask_array[j] =  rng_int_dist(rng)*2 - 1  ;
			std::cout<<mask_array[j]<<std::endl;
		}
	}
	// In this routine we get the oscillator response x_i(t), where i is the sequential node
	// i=0..24
	double time=0.0;
	double oscillation_response(double Hi){
		// define a output file to store the data

		// we calculate the no of steps needed to be performed per node
		no_steps_per_node=std::round(theta / integrate::Dt);

       		// In this loop we apply a sequence of input fields from

		//we loop over the nodes
		for (int i=1; i<=no_nodes;i++){

			//std::cout<<no_steps_per_node<<"\t"<<theta<<"\t"<<integrate::Dt<<std::endl;
			//std::cout<<stor::V0<<std::endl;
			//
			// recalculate the field
			stor::V0 = Hc + dH*Hi*mask_array[i];

			// In this loop we average over a time=theta
			for (int j=0; j<no_steps_per_node; j++){
				integrate::runge_kutta(time);
				//	std::cout<<stor::V<<"\t"<<stor::V0<<"\t"<<
                        	//	stor::x_dw<<"\t"<<time*1e9<<std::endl;
				}

			// stor the position of the domain wall into array of outputs
			s_x.push_back(stor::x_dw*1e9/150);
			outputfile << std::fixed
                       << std::setprecision(6)
					   << time*1e9 << "\t"
					   << stor::V << "\t"
					   << stor::V0 << "\t"
					   << stor::x_dw*1e9 << "\t"
                       <<std::endl;
			}

	}
	// here we define new variables for the following training process
	std::vector<double> W; // in this array we store the output weights
	const double r=0.001; // rate of learning
	double y_p=0.0; // target & output weight
	const double sigma=0.00001;
	double e_p=0.0;

	double sigmoid(double x){
	return 1/(1+exp(-x));
	}

	// in this subroutine we determine the wight matrix of the output
	// the aim is to obtain a trainer capable to classify the corresponding inputs
	// the training process is implemented using gradient method detailed in Ref. "An Introduction to Neural Networks" by Kevin Gurney, p. 90
	// Dwi= rate*(t_p-y_p)x_i, where Dwi are the weiactoghts adjustements
	double training( double lr, double la, int Nepoch, std::vector<double> &input_x, std::vector<double> &input_y){
        outputfile.open("reservoir.data");
		//initialize the mask
		mask_values();

		double bias = 0.0;
		// initialize the output weight array W
		for (int z=0; z<no_nodes; z++){
			// the weight are initialized randomly between -0.5 and 0.5
			W.push_back( rng_dist(rng) - 0.5);
		}


        int epoch = 0;
		do{
            std::cout << epoch << std::endl;
			// clear s_x
			// loop over samples
			for (int t=0; t<input_x.size(); t++){

				// delete the elements of the vector
				s_x.clear();

				// calculate the response per node
				oscillation_response(input_x[t]);
				// In this loo we calculate the activation y_p;
                // Calculate y_p = \sum_j W_j S_j
				y_p = bias;
				for(int l=0; l<no_nodes; l++){
					y_p += W[l] * s_x[l];
				}

				e_p=(input_y[t]-y_p)*(input_y[t] - y_p);
				for(int l=0; l<no_nodes; l++){
					//re-adjust the weights
					W[l] += lr*(input_y[t] - y_p)*s_x[l];
					bias += lr*(input_y[t] - y_p);
				}

				std::cout << std::fixed
                    << std::setprecision(0)
                    << t << "\t"
                    << std::setprecision(6)
                    << input_x[t] << "\t"
                    << input_y[t] << "\t"
                    << e_p << "\t"
                    << y_p << "\t"
                    << W[0] << "\t"
                    << bias
                    << std::endl;

			}
            epoch++;
            if( epoch > Nepoch) break;
		}
		while(e_p>sigma);
	//close the file
	outputfile.close();

	}

	double classification(){
		// print the weights values
		std::cout<<"Print the values of weights:"<<std::endl;
		for (int k =0; k<no_nodes;k++){
			std::cout<<W[k]<<std::endl;
		}
		std::vector <double> H_class{7.5,17.5,20,15,5,-12,-20,-12};
		// loop over test values of H
		for (int t=0; t<H_class.size(); t++){
			// delete the elements of the vector
                        s_x.clear();
			// calculate the response per node
                        oscillation_response(H_class[t]);

			//loop over the nodes and sum the x_ki
			for (int z=0; z<no_nodes; z++){
				y_p += W[z] * s_x[z];
			}

			std::cout<<y_p<<std::endl;
			if(y_p<0.5){
				std::cout<< "You got a sine input!"<<std::endl;}

		}

	return 1;
	}//end of classification function


	struct rc_inputs_t {
		std::string filename;
		double T;
		int Nv;
		int Ne;
		double lr;
		double la;
		double H0;
		double dH;
	};


	void get_input_data(std::string &filename, std::vector<double> &input_x, std::vector<double> &input_y)
	{
		std::ifstream input (filename.c_str());
		if(input.is_open() ){
			double x = 0.0, y = 0.0;
			while ( input >> x >> y ) {
				input_x.push_back(x);
				input_y.push_back(y);
			}
		}
	}


    int run()
    {
        std::vector<double> input_x;
        std::vector<double> input_y;

        // Input is read and stored in a map class
        input_map_t rc_inputs;
        rc_inputs.read_file("rc_input");

        // Print out the input values that have been recognised
        std::cout << "Stored values: " << std::endl;
        rc_inputs.print();

        Hc = rc_inputs.get<double>("H0");
        dH = rc_inputs.get<double>("dH");
        no_nodes = rc_inputs.get<int>("Nv");
        tau = rc_inputs.get<double>("T");
        theta = tau/no_nodes;

        std::cout << "tau = " << tau << std::endl;
        std::cout << "theta = " << theta << std::endl;
        std::cout << "Number of nodes = " << no_nodes << std::endl;

        // Access input values through get function with type template
        std::string filename = rc_inputs.get<std::string>("file");
        reservoir::get_input_data(filename, input_x, input_y);


        double lr = rc_inputs.get<double>("lr"); // input learning rate
        double la = rc_inputs.get<double>("la"); // input learning momentum
        int Nepoch = rc_inputs.get<double>("Ne"); // input number of epochs
        reservoir::training(lr, la, Nepoch, input_x, input_y);
        reservoir::classification();

        return 1;
    }


}
