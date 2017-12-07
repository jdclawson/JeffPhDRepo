#from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

T=0

FullData=np.genfromtxt("data_clipped.csv",delimiter=',')

Data_Series_size=FullData.shape[0]/2
Time=FullData.shape[1]

t=np.arange(Data_Series_size)

Populationgrowth=np.zeros((Data_Series_size,Time))
DebtPercent=np.copy(Populationgrowth)

for i in xrange(FullData.shape[0]):
    if i%2==0:
        Populationgrowth[i/2,:]=FullData[i,:]

    else:
        DebtPercent[i/2,:]=FullData[i,:]

#print Populationgrowth

#print DebtPercent

indextoyear={}
for i in xrange(Time):
    indextoyear[i]=2002+i

for T in xrange(Time):
    plt.ylim(0,200)
    plt.xlim(-2,3)
    plt.title ("Debt vs Population Growth in "+str(indextoyear[T]))
    plt.scatter(Populationgrowth[:,T],DebtPercent[:,T],c=t)
    plt.ylabel('Debt Percentage of GDP')
    plt.xlabel('Annual Population Growth')
    exportname=str(indextoyear[T])+".png"
    plt.savefig(exportname)
    plt.clf()
    #plt.show()
