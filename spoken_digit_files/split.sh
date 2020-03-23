#!/bin/bash
rm Training*
rm Validation*

split -l 50 -d  Processed_spoken_signal.txt 
for i in $(seq 0 1 9)
do	
	cat x0$i | head -n 25>> Training.data
	cat x0$i | tail  -n 25 >> Validation.data
done

rm x*
split -l 50 -d Y_vec.txt

for i in $(seq 0 1 9)       
do 
	cat x0$i | head -n 25 >> Trainingy.data
	cat x0$i | tail  -n 25 >> Validationy.data   
done 

