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
#include <string>

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
    array_t <1,double> MASK;

    // in this function we generate random number for the mask which can be either -1 or +1
    void mask_values(){
    	MASK.assign( 1024*no_nodes, 0.0);
        mask_array.assign( no_nodes, 0);
        for (int n=0; n<no_nodes; n++){
            mask_array[n] =  rng_int_dist(rng)*2.0 - 1.0;
            std::cout<<mask_array[n]<<std::endl;
	    for (int j=0; j<1024;j++){
		MASK(n*1024 +j )= (rng_int_dist(rng)*2.0 - 1.0);//*0.1;
	    }
        }

    }

    // In this routine we get the oscillator response x_i(t), where i is the sequential node
    // i=0..24
    void oscillation_response(double Hi){
	double time=0.0;
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


    void generate_signal_multi_dw( array_t<2, double> &mask, array_t<2,double> &Signal, std::vector<double> &input_x, const char* file)
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
            Signal(t,0) = 1.0;

            for (int i=0; i<no_nodes;i++)
            {
                //store the average position of the DW
                double average_position=0.0;

                // recalculate the field
                for( int j = 0; j < stor::Nwires; j++)
                    stor::V0_mdw[j] = Hc + dH*(0.001 + input_x[t])*mask(j,i);


                // In this loop we average over a time=theta
                for (int j=0; j<no_steps_per_node; j++){
                    integrate::multi_dw::runge_kutta(stor::x_coord, stor::phi_coord, time, integrate::Dt);
                    //average_position+= stor::x_dw*stor::x_dw*1e18;
                    //if (j%100 == 99) outstream << time*1e9 << "\t" << stor::x_dw*1e7 << "\t" << stor::V0 << std::endl;
                }

                // store the position of the domain wall into array of outputs
                //Signal(t,i) = (sqrt(average_position/no_steps_per_node));
                for( int j = 0; j < stor::Nwires; j++)
                    Signal(t,(i+1)+j*no_nodes) = (stor::x_coord[j]/1e-7);

                if ( outstream.is_open() )
                    for( int j = 0; j < stor::Nwires; j++)
                        outstream << Signal(t, i+1 + j*no_nodes) << "\t";
            }
            outstream << std::endl;
        }
        outstream.close();
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
            //outstream << t << "\t" << input_x[t] << "\t";

            for (int i=0; i<no_nodes;i++)
            {
                //store the average position of the DW
                double average_position=0.0;

                // recalculate the field
                stor::V0 = Hc + dH*(0.001 + input_x[t])*mask_array[i];

                // In this loop we average over a time=theta
                for (int j=0; j<no_steps_per_node; j++){
                    integrate::runge_kutta(time);
                    //average_position+= stor::x_dw*stor::x_dw*1e18;
                    //if (j%100 == 99) outstream << time*1e9 << "\t" << stor::x_dw*1e7 << "\t" << stor::V0 << std::endl;
                }

                // store the position of the domain wall into array of outputs
                Signal(t,i) = (sqrt(average_position/no_steps_per_node));
                // Signal(t,i) = (stor::x_dw/1e-7);
                if ( outstream.is_open() )
                    outstream << Signal(t,i) << "\t";
            }
            //outstream << std::endl;
        }
        outstream.close();
    }

    extern "C" {
        extern int dgesv_( int*, int*, double*, int*, int*, double*, int*, int*);
        extern int dgelss_( int*, int*, int*, double*, int*, double*, int*, double*, double*, int*, double*, int*, int*);
    }

    void linear_regression2( const int Nout, array_t<2,double> &Weights, array_t<2,double> &S, array_t<2,double> &input_y, double regression_factor)
    {

        int N = S.size(1);
        Weights.assign(Nout, N, 0.0);

        array_t<2,double> STS;
        STS.assign( N, N, 0.0);

        array_t<2,double> STY;
        STY.assign( N, Nout, 0.0);

        // Regularisation param
        double alpha = regression_factor;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                STS(i,j) = 0.0;
                for( int k = 0; k < S.size(0); k++) {
                    // Store the transpose
                    STS(i,j) = STS(i,j) + S(k,j) * S(k,i);
                }
            }
            // STS += alpha*I
            STS(i,i) = STS(i,i) + alpha*alpha;
        }

        std::cout<<"StS sizes are .... mxn :"<<"\t"<<STS.size(0)<<"\t"<<STS.size(1)<<std::endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Nout; j++) {
                STY(i,j) = 0.0;
                for( int k = 0; k < S.size(0); k++) {
                    STY(i,j) = STY(i,j) +  S(k,i) * input_y(k,j);
                }
            }
        }

        int NRHS = Nout;
        int LDA = N;
        int *IPIV = new int[N];
        int LDB = N;
        int INFO;

        double *SV = new double[N];
        double RCOND = 1e-10;
        int RANK;
        double *WORK = new double[6*N];
        int LWORK = 6*N;


        //dgesv_( &N, &NRHS, &STS(0,0), &LDA, IPIV, &STY(0,0), &LDB, &INFO);
        dgelss_( &N, &N, &NRHS, &STS(0,0), &LDA, &STY(0,0), &LDB, SV, &RCOND, &RANK, WORK, &LWORK, &INFO);

        // Check if linear solver has completed correctly
        if( INFO == 0) {
            // Store the result as the transposed weights
            for( int i = 0; i < Nout; i++)
                for (int j = 0; j < N; j++)
                    Weights(i,j) = STY(j,i);
        } else {
            std::cerr << "Lapack returned INFO != 0: INFO = " << INFO << std::endl;
        }

        delete [] IPIV;
        delete [] SV;
        delete [] WORK;

    }

    void linear_regression( const int Nout, array_t<2,double> &Weights, array_t<2,double> &S, array_t<2,double> &input_y, double regression_factor)
    {

        int N = S.size(1);
        Weights.assign(Nout, N, 0.0);

        array_t<2,double> STS;
        STS.assign( N, N, 0.0);

        array_t<2,double> STY;
        STY.assign( N, Nout, 0.0);

        // Regularisation param
        double alpha = regression_factor;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                STS(i,j) = 0.0;
                for( int k = 0; k < S.size(0); k++) {
                    // Store the transpose
                    STS(i,j) = STS(i,j) + S(k,j) * S(k,i);
                }
            }
            // STS += alpha*I
            STS(i,i) = STS(i,i) + alpha*alpha;
        }

        std::cout<<"StS sizes are .... mxn :"<<"\t"<<STS.size(0)<<"\t"<<STS.size(1)<<std::endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Nout; j++) {
                STY(i,j) = 0.0;
                for( int k = 0; k < S.size(0); k++) {
                    STY(i,j) = STY(i,j) +  S(k,i) * input_y(k,j);
                }
            }
        }

        int NRHS = Nout;
        int LDA = N;
        int *IPIV = new int[N];
        int LDB = N;
        int INFO;

        dgesv_( &N, &NRHS, &STS(0,0), &LDA, IPIV, &STY(0,0), &LDB, &INFO);

        // Check if linear solver has completed correctly
        if( INFO == 0) {
            // Store the result as the transposed weights
            for( int i = 0; i < Nout; i++)
                for (int j = 0; j < N; j++)
                    Weights(i,j) = STY(j,i);
        } else {
            std::cerr << "Lapack returned INFO != 0: INFO = " << INFO << std::endl;
        }

        delete [] IPIV;

    }


    void linear_regression( const int Nout, array_t<2,double> &Weights, array_t<2,double> &S, std::vector<double> &input_y)
    {

        int N = S.size(1);
        Weights.assign(Nout, N, 0.0);

        array_t<2,double> STS;
        STS.assign( N, N, 0.0);

        array_t<2,double> STY;
        STY.assign( N, Nout, 0.0);

        // Regularisation param
        double alpha = 0.001;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for( int k = 0; k < S.size(0); k++) {
                    // Store the transpose
                    STS(i,j) = STS(i,j) + S(k,j) * S(k,i);
                }
            }
            // STS += alpha*I
            STS(i,i) = STS(i,i) + alpha*alpha;
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Nout; j++) {
                for( int k = 0; k < S.size(0); k++) {
                    STY(i,j) = STY(i,j) + input_y[k] * S(k,i);
                }
            }
        }

        int NRHS = Nout;
        int LDA = N;
        int *IPIV = new int[N];
        int LDB = N;
        int INFO;

        dgesv_( &N, &NRHS, &STS(0,0), &LDA, IPIV, &STY(0,0), &LDB, &INFO);

        // Check if linear solver has completed correctly
        if( INFO == 0) {
            // Store the result as the transposed weights
            for( int i = 0; i < Nout; i++)
                for (int j = 0; j < N; j++)
                    Weights(i,j) = STY(j,i);
        } else {
            std::cerr << "Lapack returned INFO != 0: INFO = " << INFO << std::endl;
        }

        delete [] IPIV;

    }
    int accuracy_spoken(array_t<2,double> &pred, array_t<2,double> &y_vec){
        int Ncorrect=0;
        int Nfail=0;

        // check dimensionality if the arrays
        if (( pred.size(0)!=y_vec.size(0 ))&& ( pred.size(0)!=y_vec.size(0) ) )
        {
            std::cerr<<"Conflict matrix sizes"<<std::endl;

        }
        for ( int i = 0; i < pred.size(0); i++) {
            int pind=0, yind=0;
            double pmax = pred(i,0);
            double ymax = y_vec(i,0);
            for( int j = 1; j < pred.size(1); j++) {
                pind = (pmax < pred(i,j)) ? j : pind;
                pmax = (pmax < pred(i,j)) ? pred(i,j) : pmax;
                yind = (ymax < y_vec(i,j)) ? j : yind;
                ymax = (ymax < y_vec(i,j)) ? y_vec(i,j) : ymax;
            }
            std::cout << i << "  " << yind << "  " << pind << "  " << ymax << "  " << pmax << std::endl;
            Ncorrect += (yind == pind) ? 1 : 0;
        }


        return Ncorrect;

    }

    int accuracy( std::vector<double> pred, std::vector<double> &corr)
    {
        int Ncorrect = 0;
        for( int i = 0; i < corr.size(); i++)
            Ncorrect += ( ((pred[i] > 0.5) ? 1 : 0) == int(corr[i])) ? 1 : 0;

        return Ncorrect;
    }

    void linear_model_spoken( array_t<2,double> &pred, array_t<2,double> &Signal, array_t<2,double> &Weights)
    {
        //std::cout<<"W sizes are  ... m x n:"<<"\t"<<Weights.size(0)<<"\t"<<Weights.size(1)<<std::endl;
        //pred.assign(Signal.size(0), Weights.size(0), 0.0);
        for( int i = 0; i < pred.size(0); i++) {
            for ( int j = 0; j < pred.size(1); j++) {
                pred(i,j) = 0.0;
                for( int k = 0; k < Weights.size(1); k++) {
                    pred(i,j) += Weights(j,k) * Signal(i,k);
                }
            }
        }

        // transpose the matrix
        /*array_t <2,double> Weights_t;
        Weights_t.assign(Weights.size(1), Weights.size(0), 0.0);
        for(int i,j=0; i<Weights.size(1), j<Weights.size(0); i,j++){
            Weights_t(i,j)=Weights(j,i);
        }
        std::cout<<"Sizes of Weights_t :"<<"\t"<<Weights_t.size(0)<<std::endl;
        */
        //calculate the prediction matrix
        /*if(Weights_t.size(0)!=Signal.size(1)) std::cerr<<"Conflict matrix sizes!"<<std::endl;
          for (int i=0; i<Signal.size(0);i++){
          for(int j=0; j<Weights_t.size(1); j++){
          double sum=0.0;
          for( int k=0; k<Signal.size(1); k++){
          sum+=Signal(i,k)*Weights_t(k,j);}
          pred(i,j)=sum;
          }
          }*/


        //std::cout<<"Pred size mxn.................. : "<<"\t"<<pred.size(0)<<"\t"<<pred.size(1)<<std::endl;
    }

    void linear_model( std::vector<double> &pred, array_t<2,double> &Signal, array_t<2,double> &Weights)
    {
        pred.assign(Signal.size(0), 0.0);
        for( int i = 0; i < Signal.size(0); i++) {
            for( int j = 0; j < Weights.size(1); j++) {
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
        for( int i = 0; i < no_nodes; i++){
            std::cout << Weights(0,i) << std::endl;
	}
        std::vector<double> pred;
        linear_model(pred, Signal, Weights);

        int Ncorrect = accuracy(pred, input_y);

        std::cout << "Number correct = " << Ncorrect << " out of " << input_y.size() << std::endl;

        Signal.assign( valid_x.size(), stor::Nwires*no_nodes, 0.0);

        generate_signal( Signal, valid_x, "Valid_Signal.out");

        pred.clear();
        linear_model( pred, Signal, Weights);

        Ncorrect = accuracy(pred, valid_y);


        std::cout << "Validation: Number correct = " << Ncorrect << " out of " << valid_y.size() << std::endl;
        for ( int i = 0; i < valid_x.size(); i++)
            std::cout << i << "\t" << valid_x[i] << "\t" << valid_y[i] << "\t" << pred[i] << "\t" << (( ((pred[i] > 0.5) ? 1 : 0) == int(input_y[i])) ? 1 : 0) << std::endl;

        //for ( int i = 0; i < valid_y.size(); i++)
        //    std::cout << i << "\t" << valid_x[i] << "\t" << valid_y[i] << "\t" << pred[i] << std::endl;
        return Ncorrect;
    }

    double multi_dw_batch_training( std::vector<double> &input_x, std::vector<double> &input_y, std::vector<double> &valid_x, std::vector<double> &valid_y)
    {
        //initialize the mask
        array_t<2, double> mdw_mask;
        mdw_mask.assign( stor::Nwires, no_nodes, 0.0);
        std::cout << "Mask : " << std::endl;
        for( int i = 0; i < stor::Nwires; i++) {
            std::cout << "Wire " << i << " | " << std::fixed << std::setprecision(6);
            for (int j=0; j<no_nodes; j++){
                mdw_mask(i,j) =  rng_dist(rng)*2.0 - 1.0  ;
                std::cout << mdw_mask(i,j) << "\t";
            }
            std::cout << std::endl;
        }

        const int Nout = 1;
        array_t<2,double> Signal;
        array_t<2,double> Weights;

        Signal.assign( input_x.size(), 1 + stor::Nwires*no_nodes, 0.0);

        generate_signal_multi_dw( mdw_mask, Signal, input_x, "Signal_multi_dw.out");

        linear_regression( Nout, Weights, Signal, input_y );

        std::cout << "Bias = " << Weights(0,0) << std::endl;
        std::cout << "Weights = " << std::endl;
        for( int j = 0; j < stor::Nwires; j++){
            std::cout << "Wire " << j << " | " << std::fixed << std::setprecision(6);
            for( int i = 0; i < no_nodes; i++){
                std::cout << Weights(0,i+1+j*no_nodes) << "\t";
            }
            std::cout << std::endl;

        }

        std::vector<double> pred;
        linear_model(pred, Signal, Weights);

        int Ncorrect = accuracy(pred, input_y);

        std::cout << "Number correct = " << Ncorrect << " out of " << input_y.size() << std::endl;

        Signal.set_all(0.0); //assign( valid_x.size(), no_nodes, 0.0);

        generate_signal_multi_dw( mdw_mask, Signal, valid_x, "Valid_Signal_multi_dw.out");

        pred.clear();
        linear_model( pred, Signal, Weights);

        Ncorrect = accuracy(pred, valid_y);


        std::cout << "Validation: Number correct = " << Ncorrect << " out of " << valid_y.size() << std::endl;
        std::ofstream outstream;
        outstream.open("accuracy_valid.data");
        for ( int i = 0; i < valid_x.size(); i++)
          outstream << i << "\t" << valid_x[i] << "\t" << valid_y[i] << "\t" << pred[i] << "\t" << (pred[i] > 0.5 ? 1 : 0) << std::endl;
        //(( ((pred[i] > 0.5) ? 1 : 0) == int(input_y[i])) ? 1 : 0) << std::endl;
        /*
          for ( int i = 0; i < valid_y.size(); i++)
          std::cout << i << "\t" << valid_x[i] << "\t" << valid_y[i] << "\t" << pred[i] << std::endl;
          */
        outstream.close();
        return Ncorrect;
    }


    double bias;
    // the aim is to obtain a trainer capable to classify the corresponding inputs
    // the training process is implemented using gradient method detailed in Ref. "An Introduction to Neural Networks" by Kevin Gurney, p. 90
    // Dwi= rate*(t_p-y_p)x_i, where Dwi are the weiactoghts adjustements
    void training( double lr, double la, int Nepoch, std::vector<double> &input_x, std::vector<double> &input_y){
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

	void read_spectogram( array_t <2,double> &sp_sig){

		//std::ifstream file("spoken_digit_files/Nicolas.dat");
		std::ifstream file("Signal_in.data");
		if(!file) {
			std::cerr<<"Failed to open the spectogram file!"<<std::endl;
		}
		/*int num_lines=0;
		std::string line;
		while (std::getline(file, line))
			++num_lines;
		std::cout<<"No of lines is:" <<num_lines<<std::endl;

		if(num_lines==0) std::cerr<<"The spectogram is empty or is not appropriate format! "<<std::endl;
		*/
		sp_sig.assign( 2000, 1025, 1.0 );
		//loop over the rows and columns
		for (int i=0; i<2000; i++){
			for (int j=0; j<1025; j++){
				file>>sp_sig(i,j);
				//std::cout<<sp_sig(i,j)<<"\t"<<std::endl;
			}
		}
		std::cout<<"Sizes of the matrix sp_sig ........."<<sp_sig.size(0)<<"\t"<<sp_sig.size(1)<<std::endl;
		file.close();
	}

	void get_signal_digit(array_t<2,double> &sp_sig, array_t <2,double> &Xij){
		int no_steps_per_node=std::round(theta / integrate::Dt);
		Xij.assign(2000, 1024, 0.0);
		int time_steps=16;
		double time=0.0;
		std::ofstream file_proc;
		file_proc.open("Processed_spoken_signal.txt");
		if(!file_proc) std::cerr<<"Failed to open output file"<<std::endl;
		//////////////////// TODO: This part of the simulation can be paralelized ///////////////
		for (int i=0; i<Xij.size(0); i++){
			stor::x_dw=0.0, time=0.0;

			//for (int k =0; k<64; k++){
			for (int k =0; k<1024; k++){
				double avr_pos=0.0;

				//loop over the time intervals
				//for (int n=1; n<=time_steps; n++){

				//avr_pos+=sp_sig(i,k*16+n-1);

				//stor::V0 = Hc + dH*sp_sig(i,k*16+n-1); //*MASK(1024*n+k);
				stor::V0 = Hc + dH*sp_sig(i,k); //*MASK(1024*n+k);
				//std::cout<<"V0 is ................="<<stor::V0<<"\t"<<"spoken signal"<<sp_sig(i,64*(n-1)+n-1)<<std::endl;

					for (int t=0; t<no_steps_per_node; t++){

						integrate::runge_kutta(time);
						avr_pos += stor::x_dw*stor::x_dw*1e18;
					}
				//}
			//Xij(i,k) = sqrt(avr_pos / (time_steps*no_steps_per_node));
			//Xij(i,k) = (stor::x_dw*stor::x_dw/time_steps)*1e18;
			Xij(i,k) = sqrt(avr_pos)/time_steps;

			}

		}
		for (int i=0; i<Xij.size(0); i++){
			for(int j=0; j<Xij.size(1); j++){
				file_proc<<Xij(i,j)<<"\t";
			}
			file_proc<<"\n";
		}


	}

    double MSE( array_t<2,double> &pred, array_t<2,double> &desired)
    {
        double error = 0.0;
        for( int i = 0; i < pred.size(0); i++){
            for( int j = 0; j < pred.size(1); j++){
                error += 0.5*( pred(i,j) - desired(i,j))*( pred(i,j) - desired(i,j));
            }
        }
        return error/double(pred.size(0));
    }


    double spoken_training(double &regression_factor){

        mask_values();

        // in this array we store the input signal from spectogram
        array_t<2, double> sp_sig;
        read_spectogram(sp_sig);

        array_t<2,double> Xij;
        array_t<2,double> Weights;
        array_t<1,double> Y_digit;
        array_t<2,double> Y_vec;


        int Nsamples = sp_sig.size(0);

        Y_digit.assign(Nsamples, 0.0);
        Y_vec.assign(Nsamples, 10, 0.0);

        for (int i = 0; i < Nsamples; i++){
            Y_digit(i) = sp_sig(i,1024);
            Y_vec(i, Y_digit(i)) = 1.0;
        }

        /*for (int i = 0; i < Nsamples; i++){
            for (int j=0; j<10; j++){
                std::cout <<Y_vec(i,j)<<"\t";
            }
            std::cout<<std::endl;
        }*/

        std::cout << "Nsamples = " << Nsamples << std::endl;
       /* Xij.assign( 500, 1024*no_nodes + 1, 0.0);
        for (int i = 0; i < Nsamples; i++){
            Xij(i, 1024*no_nodes) = 1.0;
            for ( int j = 0; j < 1024*no_nodes; j++){
                Xij(i,j) = sp_sig(i,j);
                //std::cerr << i << "  " << j << "  " << Xij(i,j) << std::endl;
            }
            //std::cerr << std::endl;
            //std::cerr << std::endl;
        }*/
        //std::cout << "Calculating signal now......" << std::flush;
        get_signal_digit(sp_sig, Xij);

        std::cout << " Done" << std::endl;
        std::cout << "Starting training now ....." << std::flush;

        //Weights.assign(10, (sp_sig.size(1)-1)*no_nodes, 0.0);
        linear_regression2( 10, Weights, Xij, Y_vec, regression_factor);

        /*
        Weights.assign( 10, Xij.size(1), 0.0);
        for (int i=0; i<Weights.size(0);i++){
            for (int j=0; j<Weights.size(1);j++){
                Weights(i,j) = 0.1*(rng_dist(rng)*2.0 - 1.0);
                //std::cout<<Weights(i,j)<<"\t";
            }
            //std::cout<<"\n";
        } */

        array_t <2,double> Pred_vec;

        Pred_vec.assign(Nsamples, 10, 0.0);


        linear_model_spoken(Pred_vec, Xij, Weights);

        double err = 0.0;
        err = MSE( Pred_vec, Xij);

        std::cout << "Mean Squared Error = " << err << std::endl;

        for (int i=0; i<Pred_vec.size(0);i++){
            for (int j=0; j<Pred_vec.size(1);j++){
                std::cout << std::fixed << std::setprecision(4) << Pred_vec(i,j) << "  ";
            }
            std::cout<<"\n";
        }

        //for (int i= 0; i < Weights.size(0); i++){
        //       for(int j = 0; j<Weights.size(1); j++){
        //   std::cout << i << "  " << j << "  " << Weights(i,j) << std::endl;
        //}}

        int Ncorrect = accuracy_spoken(Pred_vec, Y_vec);
        std::cout<< "No of correct predictions:  " << Ncorrect << " out of " << Pred_vec.size(0) << std::endl;

    }

    int run_spoken_recognition(){

        input_map_t rc_inputs;
        rc_inputs.read_file("rc_input");
        std::cout << "Stored values: " << std::endl;
        rc_inputs.print();
        int seed = rc_inputs.get<int>("seed");
        rng.seed(seed);
        Hc = rc_inputs.get<double>("H0");
        dH = rc_inputs.get<double>("dH");
        no_nodes = rc_inputs.get<int>("Nv");
        tau = rc_inputs.get<double>("T");
        theta = tau/no_nodes;
        double regression_factor = rc_inputs.get<double>("lambda");

        std::cout << "tau = " << tau << std::endl;
        std::cout << "theta = " << theta << std::endl;
        std::cout << "Number of nodes = " << no_nodes << std::endl;
        std::cout << "Regression factor, lambda = "<<regression_factor << std::endl;

        spoken_training(regression_factor);
        return 0;
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

        /*
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
        */
        multi_dw_batch_training( x_train, y_train, x_valid, y_valid);

        //reservoir::training(lr, la, Nepoch, x_train, y_train);
        //reservoir::classification(x_valid, y_valid);

        return 1;
    }




    int run_field_sequence()
    {

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

        int Nv = rc_inputs.get<int>("Nv");
        int Ns = rc_inputs.get<int>("Ns");

        array_t<2,double> input_x;
        array_t<2,double> Signal;

        input_x.assign(Ns, Nv, 0.0);
        Signal.assign(Ns, Nv, 0.0);

        std::ifstream input (filename.c_str());
        if(input.is_open() ){
            for (int i = 0; i < Ns; i++) {
                for ( int j = 0; j < Nv; j++) {
                    double x =0;
                    input >> x;
                    input_x(i,j) = x;
                }
            }
            input.close();
        }



        std::string outfilename = rc_inputs.get<std::string>("outfile");
        std::ofstream outstream;
        outstream.open(outfilename.c_str());

        // we calculate the no of steps needed to be performed per node
        no_steps_per_node=std::round(theta / integrate::Dt);
        std::cout << "Steps per node = " << no_steps_per_node << std::endl;

        double time = 0.0;
        for (int t=0; t<Ns; t++)
        {
            //outstream << t << "\t" << input_x[t] << "\t";

            for (int i=0; i<no_nodes;i++)
            {
                //store the average position of the DW
                double average_position=0.0;

                // recalculate the field
                stor::V0 = Hc + dH*input_x(t,i);

                // In this loop we average over a time=theta
                for (int j=0; j<no_steps_per_node; j++){
                    integrate::runge_kutta(time);
                    average_position+= stor::x_dw*stor::x_dw*1e18;
                    //if (j%100 == 99) outstream << time*1e9 << "\t" << stor::x_dw*1e7 << "\t" << stor::V0 << std::endl;
                }

                // store the position of the domain wall into array of outputs
                Signal(t,i) = sqrt(average_position/no_steps_per_node);
                // Signal(t,i) = (stor::x_dw/1e-7);
                if ( outstream.is_open() )
                    outstream << Signal(t,i) << "\t";
            }
            outstream << std::endl;
        }
        outstream.close();

        return 1;
    }

    int run_transient()
    {

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

        int Nv = no_nodes;
        int Nd = rc_inputs.get<int>("Nd");
        int Ns = rc_inputs.get<int>("Ns");

        array_t<2,double> input_x;
        array_t<2,double> Signal;

        input_x.assign(Ns, Nd, 0.0);
        Signal.assign(Ns, Nd*Nv, 0.0);

        std::ifstream input (filename.c_str());
        if(input.is_open() ){
            for (int i = 0; i < Ns; i++) {
                for ( int j = 0; j < Nd; j++) {
                    double x =0;
                    input >> x;
                    input_x(i,j) = x;
                }
            }
            input.close();
        }



        std::string outfilename = rc_inputs.get<std::string>("outfile");
        std::ofstream outstream;
        outstream.open(outfilename.c_str());

        // we calculate the no of steps needed to be performed per node
        no_steps_per_node=std::round(theta / integrate::Dt);
        std::cout << "Steps per node = " << no_steps_per_node << std::endl;

        double time = 0.0;
        for (int t=0; t<Ns; t++)
        {
            //outstream << t << "\t" << input_x[t] << "\t";

            // loop over the input dimensions
            for (int i=0; i<Nd;i++)
            {
                // recalculate the field
                stor::V0 = Hc + dH*input_x(t,i);
		stor::x_dw = 0.0;
		stor::phi_dw = 0.0;

                for ( int n=0; n < Nv; n++){

                    // In this loop we average over a time=theta
                    for (int j=0; j<no_steps_per_node; j++){
                        integrate::runge_kutta(time);
                        //if (j%100 == 99) outstream << time*1e9 << "\t" << stor::x_dw*1e7 << "\t" << stor::V0 << std::endl;
                    }

                    // store the position of the domain wall into array of outputs
                     Signal(t,Nv*i + n) = (stor::x_dw/1e-7);
                    if ( outstream.is_open() )
                        outstream << Signal(t,Nv*i + n) << "\t";
                }

            }
            outstream << std::endl;
        }
        outstream.close();

        return 1;
    }
}
