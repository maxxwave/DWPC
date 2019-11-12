# DWPC
This code is dedicated to simulate a 1D domain-wall propagation into a nano-strip with two notches or anti-notches wich give power law pinning energy. The domain wall is driven by an oscillating magnetic field between the notches. The roles of notches/anti-notches is to give a non-linear oscillation response similar to a Duffin oscillator. 

As an extension of this, the dynamics of the DW can be applied to reproduce the multiplexing of a signal into a Recurrent Neural Netword (RNN) in order to obtain a classification of different inputs signals.

This code is developed by Razvan Ababei and Matt Ellis in collaboration with Tom Hayward at University of Sheffield, UK.


# Installing DWPC
DWPC requires to be run on Linux OS using g++ compiler. The installation can be done by cloning github repository as follow:

"git clone https://github.com/maxxwave/DWPC.git"

Once this is downloaded in your local Linux computer you can run the installation script following:

"cd DWPC"
&
"sh makefile.sh"

# Using DWPC
The parameters can be changed from the input file. An example comes with the package. The user must be aware of respecting the format of the input file. No space is allowed between the name of the variable and "=" and between "=" and the value.

Note: All variables are set in SI units

