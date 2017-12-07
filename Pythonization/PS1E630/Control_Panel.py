from __future__ import division #Forces regular division for integers
import time #Package for timing things
import numpy as np  #Matrix algebra package
import Econ630Algorithms as AUX #Algorithms Class
import Functions as F #Auxiliary Functions
from scipy.interpolate import interp1d

'''
This is a python version of code prepared by various class members of Econ 630 at
the University of Maryland. It's designed to be as readable and straightfoward as possible and 
give the first time learner an example to go off of. I'm going to assume that you don't have
much experience reading Python code, so I'll comment on specific differences that you may or
may not be aware of

Python Conversion (and some original MATLAB code) By: Jeff Clawson


'''

#Parameters
A = 1 #Original productivity level
alpha = .27 #Capital share of production
beta = .994 #Paitence Parameter
delta = .011 #Depreciation rate
eta = 2 #Coefficient of relative risk aversion
#Theta is set so that h_ss is .33
N = 15 #Policy Function iteration update parameter
n= 300 #Number of grid points for capital
multi_gridpts = [50,150,300]
multi_gridptsint = [50,300]
A_Tran = [1,2]
T = 500 #Number of periods of shooting algorithm
tol = 1e-12 #Tolerance for convergence
shoottol=1e-7
upper = 1.2 #Factors to scale up the upper bound for the capital grid
lower = .5 #Factors to scale down the lower bound for the capital grid
maxIter = 3000 #Maximum number of iterations

params=(A, alpha, beta, delta, eta,N,tol,maxIter,upper,lower,shoottol)


#Levers, True if you want it to run, false otherwise
#VFI
NaiveVFI = False
NormalizeSSVFI = False 
Multi_Grid = False
MGMonotonicity = False
MGConcaveBS = False
MGCBSMcQueenPort = True
LinInterpol = True

#Policy Function Iteration
PFIAccelerator = False
PFIPoormans = True

#Shooting Algorithm
ShootingSS = False

#Progress
showprogress = True
showgraphs = True

#Initialize the Algorithm Object
Model=AUX.Econ_Algorithm(params)

Euler_Params = (A, alpha, beta, delta, eta,Model.theta, n)

if NaiveVFI:
    start=time.time()
    nVFI_guess = np.zeros(n)
    Model.VFI(showprogress,n,nVFI_guess)
    NaiveVFItime=time.time()-start

    F.plotVarPol(Model.nVFIkgrid,Model.nVFIValue,Model.nVFIPolicy)

    F.CalcResid(Model.nVFIkgrid,Model.nVFIValue,Model.nVFIIndex,Model.nVFIHmat,Euler_Params)

if NormalizeSSVFI:

    NSSVFI_guess = np.ones(n)*Model.kss
    if NaiveVFI:
        Model.VFI(showprogress,n,NSSVFI_guess,scale=(1-beta),reuse=Model.nVFIHmat)
    else:
        Model.VFI(showprogress,n,NSSVFI_guess,scale=(1-beta))

    F.plotVarPol(Model.NSSVFIkgrid,Model.NSSVFIValue,Model.NSSVFIPolicy)

    F.CalcResid(Model.NSSVFIkgrid,Model.NSSVFIValue,Model.NSSVFIIndex,Model.nVFIHmat,\
            Euler_Params)

if Multi_Grid:



    MG_guess=np.ones(multi_gridpts[0])*Model.kss

    for i in xrange(len(multi_gridpts)):
        Model.VFI(showprogress,multi_gridpts[i],MG_guess,scale=(1-beta),MG=True)

        if i!=len(multi_gridpts)-1:
            x_i=np.arange(multi_gridpts[i])
            x_i1=np.arange(multi_gridpts[i+1])

            New_Guess=np.interp(x_i1,x_i,Model.MGVFIvalue)

    Euler_Params = (A, alpha, beta, delta, eta,Model.theta, multi_gridpts[i])

    F.plotVarPol(Model.MGVFIkgrid,Model.MGVFIvalue,Model.MGVFIPolicy)

    F.CalcResid(Model.MGVFIkgrid,Model.MGVFIvalue,Model.MGVFIIndex,Model.MGVFIHmat\
            ,Euler_Params)

if MGMonotonicity:

    MGMonoguess=np.ones(multi_gridptsint[0])*Model.kss

    for i in xrange(len(multi_gridptsint)):
        Model.VFI(showprogress,multi_gridptsint[i],MGMonoguess,MG=True,Mono=True)


        if i != len(multi_gridptsint)-1:
            x_i=np.arange(multi_gridptsint[i])
            x_i1=np.arange(multi_gridptsint[i+1])

            New_Guess=np.interp(x_i1,x_i,Model.MGVFIvalue)

    Euler_Params = (A, alpha, beta, delta, eta,Model.theta, multi_gridptsint[i])
       
    F.plotVarPol(Model.MGVFIkgrid,Model.MGVFIvalue,Model.MGVFIPolicy)
    
    F.CalcResid(Model.MGVFIkgrid,Model.MGVFIvalue,Model.MGVFIIndex,Model.MGVFIHmat\
            ,Euler_Params)

if MGConcaveBS:

    MGConcaveGuess= np.ones(multi_gridptsint[0])*Model.kss

    for i in xrange(len(multi_gridptsint)):
        Model.VFI(showprogress,multi_gridpts[i],MGConcaveGuess,scale=(1-beta),\
                MG=True,Concave=True)

        if i != len(multi_gridptsint)-1:
            x_i=np.arange(multi_gridptsint[i])
            x_i1=np.arange(multi_gridptsint[i+1])

            New_Guess=np.interp(x_i1,x_i,Model.MGVFIvalue)

    F.plotVarPol(Model.MGVFIkgrid,Model.MGVFIvalue,Model.MGVFIPolicy)
    
    F.CalcResid(Model.MGVFIkgrid,Model.MGVFIvalue,Model.MGVFIIndex,Model.MGVFIHmat\
            ,Euler_Params)

if PFIAccelerator:
    PFI_guess= np.ones(n)*Model.kss


    Model.PFI(PFI_guess,n,showprogress)


if PFIPoormans:

    PFI_guess= np.ones(n)*Model.kss

    
    Model.PFI(PFI_guess,n,showprogress,poor_n=N,poormans=True)

if ShootingSS:

    Model.ShootingAlgorithm(A_Tran,T,showprogress)

    F.PlotTransPath(T,Model.cshootpath,Model.kshootpath)



        




