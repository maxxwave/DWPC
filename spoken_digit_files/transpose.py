import sys
import numpy as np
#f=open('Processed_spoken_signal.txt','r')
#f=open('Theo.dat','r')

f=open('histogram.dat','r')

xij=[]
xij=[line.split() for line in f]
xij=np.array(xij).astype('float')
xij= xij.transpose()
for i in range(len(xij)):
    print(' '.join(map(str, xij[i]))) 
f.close()
