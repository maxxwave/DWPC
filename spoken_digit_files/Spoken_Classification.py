#!/bin/python

import numpy as np 
from numpy.linalg import inv
from numpy.linalg import det

import csv
from numpy.linalg import inv
from numpy.linalg import det

f = open('Training.data', 'r') 
g = open('Trainingy.data', 'r')

h = open('Validation.data', 'r')
k = open('Validationy.data', 'r')



no_nodes=2
Nin = 1024*no_nodes
Ndigit = 10 
Nsample = 250

xij = [Nin]
xij = [line.split() for line in f]
xij = np.array(xij).astype('float')
xij = xij.transpose()
xijT = xij.transpose()

ytarget = [Ndigit]
ytarget = [line.split() for line in g]
ytarget = np.array(ytarget).astype('float')

validationij = [Nin] 
validationij = [line.split() for line in h]
validationij = np.array(validationij).astype('float')
validationij = validationij.transpose()

validationy = [Ndigit]
validationy = [line.split() for line in k] 
validationy = np.array(validationy).astype('float')

alpha=0.00005
print(ytarget)

Wout = np.matmul(xij,xijT)
Wout = np.add(Wout, alpha*np.identity(len(Wout)))
Wout = np.linalg.pinv(Wout)
Wout = np.matmul(xijT, Wout)
 
Wout = np.matmul(ytarget.transpose(),Wout)


#print("Wout sizes: ..", Wout.shape)
#print (Wout)

Ypred=np.matmul(Wout, validationij)

print (Ypred)

Ypred = Ypred.transpose()

print("Ypred sizes, Y_vec sizes: ", Ypred.shape, ytarget.shape)

def accuracy(Y_pred,Yvec):
    count = 0
    for i in range(0,250,1):
        a=Y_pred[i,:]
        b=Yvec[i,:]
        ind_label = np.where(a==a.max())
        ind_pred  = np.where(b==b.max()) 
        if ind_label==ind_pred:
            count +=1
        print(ind_label, ind_pred)
    print("No of predictions", count)

accuracy(Ypred,validationy)

