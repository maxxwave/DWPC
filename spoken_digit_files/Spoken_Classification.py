#!/bin/python

import numpy as np 
from numpy.linalg import inv
from numpy.linalg import det

import csv
from numpy.linalg import inv
from numpy.linalg import det

f = open('Processed_signal_spoken.txt', 'r') 
g = open('Y_vec.txt', 'r')


no_nodes=2
Nin = 1024*no_nodes
Ndigit = 10 
Nsample = 500


xij = [Nin]
xij = [line.split() for line in f]
xij = np.array(xij).astype('float')
xij = xij.transpose()
xijT = xij.transpose()

ytarget = [Ndigit]
ytarget = [line.split() for line in g]
ytarget = np.array(ytarget).astype('float')
alpha=0.000005
print(ytarget)



Wout = np.matmul(xij,xijT)
Wout = np.add(Wout, alpha*np.identity(len(Wout)))
Wout = np.linalg.pinv(Wout)
Wout = np.matmul(xijT, Wout)
 
Wout = np.matmul(ytarget.transpose(),Wout)


#print("Wout sizes: ..", Wout.shape)
#print (Wout)

Ypred=np.matmul(Wout, xij)

print (Ypred)

Ypred = Ypred.transpose()

print("Ypred sizes, Y_vec sizes: ", Ypred.shape, ytarget.shape)

def accuracy(Y_pred,Yvec):
    count = 0
    for i in range(0,500,1):
        a=Y_pred[i,:]
        b=Yvec[i,:]
        ind_label = np.where(a==a.max())
        ind_pred  = np.where(b==b.max()) 
        if ind_label==ind_pred:
            count +=1
        print(ind_label, ind_pred)
    print("No of predictions", count)

accuracy(Ypred,ytarget)

