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
x_init =.001
p_gridpts = 20
s_gridpts = 100

arg_in = (k_init,x_init)
vfi_in = None

#Levers
ShowProgress=True

params = (g,rho,beta,alpha,gamma,delta,sig_e,p_gridpts,disp,tol,s_gridpts)

params2 = np.array(params)
paramsdf = pd.DataFrame(params2)

print paramsdf.to_latex()

Basic = AUX.Project_Model(params)

Basic.modelsolve(arg_in,vfi_in,ShowProgress)
