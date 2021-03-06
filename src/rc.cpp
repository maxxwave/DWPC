// This file is dedicated to simulate the reservoir computing using a nano-strip DW oscillator
// Author: Razvan Ababei & Matt Ellis
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

#include "../hdr/storage.h"
#include "../hdr/calculate.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/rc.h"
#include "../hdr/input_map.h"

#include "../hdr/arrays.h"

namespace reservoir{

    	std::default_random_engine rng;
    	std::uniform_real_distribution<double> rng_dist(0.0,1.0);
    	std::uniform_int_distribution<int> rng_int_dist(0,1);

	std::ofstream outputfile;

	// Define a vector to store the amplitudes of the field
    // Hc defines the centre field amplitude
    double Hc = 1000;
    // dH defines the width of the field amplitude
    double dH = 100;

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
			mask_array[j] =  rng_int_dist(rng)*2.0 - 1.0  ;
			std::cout<<mask_array[j]<<std::endl;
		}
	}

	// In this routine we get the oscillator response x_i(t), where i is the sequential node
	// i=0..24
	double time=0.0;
	double oscillation_response(double Hi){

		// we calculate the no of steps needed to be performed per node
		no_steps_per_node=std::round(theta / integrate::Dt);

       		// In this loop we apply a sequence of input fields from

		//we loop over the nodes
		for (int i=0; i<no_nodes;i++){
			//store the average position of the DW
            double average_position=0.0;

			// recalculate the field
			stor::V0 = Hc + dH*(0.001 +Hi)*mask_array[i];

			// In this loop we average over a time=theta
			for (int j=0; j<no_steps_per_node; j++){
				integrate::runge_kutta(time);
				average_position+=stor::x_dw*stor::x_dw*1e18;


			// stor the position of the domain wall into array of outputs
			//s_x.push_back(sqrt(average_position/no_steps_per_node));
            s_x.push_back(stor::x_dw/1e-7);
			}
			outputfile << std::fixed
                       << std::setprecision(6)
					   << time*1e9 << "\t"
					   << stor::V << "\t"
					   << stor::V0 << "\t"
					   << stor::x_dw*1e9 << "\t"
                       << sqrt(average_position/(no_steps_per_node)) << "\t"
                       << Hi
                       <<std::endl;

        }
	}

	// here we define new variables for the following training process
	std::vector<double> W; // in this array we store the output weights
	const double r=0.001; // rate of learning
	double y_p=0.0; // target & output weight
	const double sigma=0.001;
	double e_p=0.0;

	double sigmoid(double x){
	    return 1/(1+exp(-x));
	}


    void generate_signal( array_t<2,double> &Signal, std::vector<double> &input_x, const char* file)
    {

        std::ofstream outstream;
        outstream.open(file);

        // we calculate the no of steps needed to be performed per node
        no_steps_per_node=std::round(theta / integrate::Dt);
        std::cout << "Steps per node = " << no_steps_per_node << std::endl;

        double time = 0.0;
        for (int t=0; t<input_x.size(); t++)
        {
            outstream << t << "\t" << input_x[t] << "\t";

            for (int i=0; i<no_nodes;i++)
            {
                //store the average position of the DW
                double average_position=0.0;

                // recalculate the field
                stor::V0 = Hc + dH*(0.001 + input_x[t])*mask_array[i];

                // In this loop we average over a time=theta
                for (int j=0; j<no_steps_per_node; j++){
                    integrate::runge_kutta(time);
                    average_position+= stor::x_dw*stor::x_dw*1e18;
                    //if (j%100 == 99) outstream << time*1e9 << "\t" << stor::x_dw*1e7 << "\t" << stor::V0 << std::endl;
                }

                // store the position of the domain wall into array of outputs
                 Signal(t,i) = (sqrt(average_position/no_steps_per_node));
                // Signal(t,i) = (stor::x_dw/1e-7);
                if ( outstream.is_open() )
                    outstream << Signal(t,i) << "\t";
            }
            outstream << std::endl;
        }
        outstream.close();
    }

    extern "C" {
        extern int dgesv_( int*, int*, double*, int*, int*, double*, int*, int*);
    }


    void linear_regression( const int Nout, array_t<2,double> &Weights, array_t<2,double> &S, std::vector<double> &input_y)
    {
        Weights.assign(Nout, no_nodes, 0.0);

        array_t<2,double> STS;
        STS.assign( no_nodes, no_nodes, 0.0);

        array_t<2,double> STY;
        STY.assign( no_nodes, Nout, 0.0);

        // Regularisation param
        double alpha = 0.1;

        for (int i = 0; i < no_nodes; i++) {
            for (int j = 0; j < no_nodes; j++) {
                for( int k = 0; k < S.size(0); k++) {
                    // Store the transpose
                    STS(i,j) = STS(i,j) + S(k,j) * S(k,i);
                }
            }
            // STS += alpha*I
            STS(i,i) = STS(i,i) + alpha*alpha;
        }

        for (int i = 0; i < no_nodes; i++) {
            for (int j = 0; j < Nout; j++) {
                for( int k = 0; k < S.size(0); k++) {
                    STY(i,j) = STY(i,j) + input_y[k] * S(k,i);
                }
            }
        }

        int N = no_nodes;
        int NRHS = Nout;
        int LDA = no_nodes;
        int *IPIV = new int[N];
        int LDB = no_nodes;
        int INFO;

        dgesv_( &N, &NRHS, &STS(0,0), &LDA, IPIV, &STY(0,0), &LDB, &INFO);

        // Check if linear solver has completed correctly
        if( INFO == 0) {
            // Store the result as the transposed weights
            for( int i = 0; i < Nout; i++)
                for (int j = 0; j < no_nodes; j++)
                    Weights(i,j) = STY(j,i);
        } else {
            std::cerr << "Lapack returned INFO != 0: INFO = " << INFO << std::endl;
        }

        delete [] IPIV;

    }

    int accuracy( std::vector<double> pred, std::vector<double> &corr)
    {
        int Ncorrect = 0;
        for( int i = 0; i < corr.size(); i++)
            Ncorrect += ( ((pred[i] > 0.5) ? 1 : 0) == int(corr[i])) ? 1 : 0;

        return Ncorrect;
    }


    void linear_model( std::vector<double> &pred, array_t<2,double> &Signal, array_t<2,double> &Weights)
    {
        pred.assign(Signal.size(0), 0.0);
        for( int i = 0; i < Signal.size(0); i++) {
            for( int j = 0; j < no_nodes; j++) {
                pred[i] += Weights(0,j) * Signal(i,j);
            }
        }
    }

    double batch_training( std::vector<double> &input_x, std::vector<double> &input_y, std::vector<double> &valid_x, std::vector<double> &valid_y)
    {
		//initialize the mask
		mask_values();

        const int Nout = 1;
        array_t<2,double> Signal;
        array_t<2,double> Weights;

        Signal.assign( input_x.size(), no_nodes, 0.0);

        generate_signal( Signal, input_x, "Signal.out");

        linear_regression( Nout, Weights, Signal, input_y);

        std::cout << "Weights = " << std::endl;
        for( int i = 0; i < no_nodes; i++)
            std::cout << Weights(0,i) << std::endl;

        std::vector<double> pred;
        linear_model(pred, Signal, Weights);

        int Ncorrect = accuracy(pred, input_y);

        std::cout << "Number correct = " << Ncorrect << " out of " << input_y.size() << std::endl;

        Signal.assign( valid_x.size(), no_nodes, 0.0);

        generate_signal( Signal, valid_x, "Valid_Signal.out");

        pred.clear();
        linear_model( pred, Signal, Weights);

        Ncorrect = accuracy(pred, valid_y);


        std::cout << "Validation: Number correct = " << Ncorrect << " out of " << valid_y.size() << std::endl;
        /*for ( int i = 0; i < valid_x.size(); i++)
            std::cout << i << "\t" << valid_x[i] << "\t" << valid_y[i] << "\t" << pred[i] << "\t" << (( ((pred[i] > 0.5) ? 1 : 0) == int(input_y[i])) ? 1 : 0) << std::endl;

        for ( int i = 0; i < valid_y.size(); i++)
            std::cout << i << "\t" << valid_x[i] << "\t" << valid_y[i] << "\t" << pred[i] << std::endl;
          */
        return Ncorrect;
    }


	double bias;
	// the aim is to obtain a trainer capable to classify the corresponding inputs
	// the training process is implemented using gradient method detailed in Ref. "An Introduction to Neural Networks" by Kevin Gurney, p. 90
	// Dwi= rate*(t_p-y_p)x_i, where Dwi are the weiactoghts adjustements
	double training( double lr, double la, int Nepoch, std::vector<double> &input_x, std::vector<double> &input_y){
        outputfile.open("reservoir.data");
		//initialize the mask
		mask_values();

		bias = 0.0;
		// initialize the output weight array W
		for (int z=0; z<no_nodes; z++){
			// the weight are initialized randomly between -0.5 and 0.5
			W.push_back( rng_dist(rng) - 0.5);
		}


        int epoch = 0;
        do{
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
                y_p = sigmoid(y_p);

                e_p=(input_y[t]-y_p)*(input_y[t] - y_p);
                for(int l=0; l<no_nodes; l++){
                    //re-adjust the weights
                    W[l] += lr*(input_y[t] - y_p)*s_x[l]*(y_p*(1-y_p));
                }
                bias += lr*(input_y[t] - y_p)*(y_p*(1.0 - y_p));

                std::cout << std::fixed
                    << std::setprecision(0)
                    << t << "\t"
                    << std::setprecision(6)
                    << input_x[t] << "\t"
                    << input_y[t] << "\t"
                    << e_p << "\t"
                    << y_p << "\t"
                    << s_x[0] << "\t"
                    << s_x[1]
                    << std::endl;

            }
            epoch++;
            if( epoch > Nepoch) break;
        }
        while(e_p>sigma);
	//close the file
	outputfile.close();

	}

	double classification(std::vector<double> &input_x, std::vector<double> &input_y){
	// Define a variable to store the average error
    double average_e=0.0;
	int Rate_success_sine=0, Rate_success_square=0;
	int count=0;
        // print the weights values
        std::cout<<"Print the values of weights:"<<std::endl;
        for (int k =0; k<no_nodes;k++){
            std::cout<<W[k]<<std::endl;
        }
       	 std::cout<<"Print the values of yk:"<<std::endl;

        // loop over test values of H
        for (int t=0; t<input_x.size(); t++){
            // delete the elements of the vector
            s_x.clear();
            // calculate the response per node
            oscillation_response(input_x[t]);
            y_p=bias;
            //loop over the nodes and sum the x_ki
            for (int z=0; z<no_nodes; z++){
                y_p += W[z] * s_x[z];
            }
            y_p = sigmoid(y_p);

            // if predict sine
            if( y_p <0.5 && input_y[t]==0.0){
	    	    Rate_success_sine++;
	        }

	        if(y_p>=0.5 && input_y[t]==1.0){
	    	    Rate_success_square++;
	        }

        count ++;

        std::cout<<"Validation Error:" << "\t"<<fabs(y_p - input_x[t])<<"\t"<<input_y[t]<<"\t"<<y_p <<"\t"<<bias <<std::endl;
        }
        std::cout<<"Rate of success" << "\t"<<(Rate_success_sine+Rate_success_square)/(input_x.size()*0.5)<<"\t"
            <<"No of sines & sqaures detected"<<  "\t"<<Rate_success_sine<<"\t"<<Rate_success_square<<std::endl;

        return 1;
    	}//end of classification function


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

        //Hc = rc_inputs.get<double>("H0");
        //dH = rc_inputs.get<double>("dH");
        int seed = rc_inputs.get<int>("seed");
        rng.seed(seed);


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

        std::vector<double> x_train, y_train, x_valid, y_valid;
        for( int i = 0; i < input_x.size()/2; i++) {
            x_train.push_back(input_x[i]);
            y_train.push_back(input_y[i]);
            x_valid.push_back(input_x[i + input_x.size()/2]);
            y_valid.push_back(input_y[i + input_x.size()/2]);
        }



        double lr = rc_inputs.get<double>("lr"); // input learning rate
        double la = rc_inputs.get<double>("la"); // input learning momentum
        int Nepoch = rc_inputs.get<double>("Ne"); // input number of epochs

        // In this section we calculate the success rate for a range of Hc and dH
        std::ofstream field_map_data;
        field_map_data.open("field_map.data");
        field_map_data<<"#"<<"HC"<<"\t"<<"dH"<<"\t" <<"N correct"<<std::endl;
        double rate_succ=0.0;
        for (int i=0; i<=50; i++)
        {
            Hc=i*50+100;
            for (int j=1; j<=i;j++)
            {
                dH=50*j;
                rate_succ=batch_training(x_train, y_train, x_valid, y_valid);
                field_map_data<< Hc<<"\t"<<dH<<"\t"<<rate_succ<<"\t"<<std::endl;
            }
        }
        field_map_data.close();

        //batch_training( x_train, y_train, x_valid, y_valid);

        //reservoir::training(lr, la, Nepoch, x_train, y_train);
        //reservoir::classification(x_valid, y_valid);

        return 1;
    }


}
