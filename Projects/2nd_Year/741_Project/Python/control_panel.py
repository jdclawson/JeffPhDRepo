import numpy as np
import modelhub as AUX
import pandas as pd


#Parameters
g = .03
gamma = 3
rho = .8
beta = .95
alpha = .35
delta = .04
sig_e = .03
disp = 3
tol = 1e-8

#Intial Guesses
k_init = .5 #Initial guess for steady state capital ratio
x_init = .0001
p_gridpts = 20
s_gridpts = 100

arg_in = (k_init,x_init)
vfi_in = None

#Levers
ShowProgress=True

params = (g,rho,beta,alpha,gamma,delta,sig_e,p_gridpts,disp,tol,s_gridpts)

Basic = AUX.Project_Model(params)

Basic.modelsolve(arg_in,vfi_in,ShowProgress)

verify=np.array([Basic.k_ss,Basic.x_ss])

