from __future__ import division
import csv
import time
import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.interpolate as interpol
import decimal
import tauchen as tn
from matplotlib import pyplot as plt
from pprint import pprint
from mpl_toolkits.mplot3d import Axes3D
from pylab import savefig


class Project_Model(object):

    def __init__(self,parameters):

        (self.g,self.rho,self.beta,self.alpha,self.gamma,self.delta,self.sig_e,\
                self.p_gridpts,self.disp,self.tol,self.s_gridpts)=parameters

    def marg_util(self,c):
        return c**-self.gamma

    def util(self,c):
        return c**(1-self.gamma)/(1-self.gamma)


    def steadystate(self,arg_in):
        k_in=arg_in[0]
        x_in=arg_in[1]

        r = self.alpha*k_in**(self.alpha-1)-self.delta
        w = (1-self.alpha)*k_in**self.alpha
        cy = w - x_in - k_in
        co = (1+self.g)*x_in+(1+r)*k_in

        resid=np.zeros(2)
        if cy<0 or co <0:
            print "Punishing fsolve for negative consumption"
            resid[0]=1e10
            resid[1]=1e10
        else:
            resid[0] = self.marg_util(cy)-self.beta*(1+r)*self.marg_util(co)
            resid[1] = k_in**self.alpha-co-cy-k_in-x_in

        #print resid

        return resid

    def get_ssvalues(self):
        r = self.alpha*self.k_ss**(self.alpha-1)
        w = (1-self.alpha)*self.k_ss**self.alpha
        cy = w - self.x_ss - self.k_ss
        co = (1+self.g)*self.x_ss+(1+r)*self.k_ss

        return r, w, cy, co

    def VFI(self,init_guess):

        dist=1

        p_grid, trans_k, dist = tn.tauchen(self.rho,self.sig_e,self.p_gridpts,self.disp)

        k_grid = np.linspace(.5*self.k_ss,1.2*self.k_ss,num=self.s_gridpts)

        #Value = young utility + beta * v_old

        v_possibilities = np.zeros((self.s_gridpts,self.e_gridpts))
        cy = np.copy(v_possibilities)

        while dist>self.tol:
            for i in xrange(self.s_gridpts):
                cy[i,:]=(1-self.alpha)*k_grid[i]**self.alpha-x-k_grid[i]

                cy[i,:][cy[i,:]<0]=-np.inf
                ucy[i,:]=self.util(cy)
                ucy[i,:][ucy[i,:]==np.inf]=-np.inf

            for i in xrange(self.s_gridpts):
                pass



        opt_value = np.zeros(self.p_gridpts)

        return opt_value


    def modelsolve(self,initss,initvfi,showprogress):

        if showprogress:
            print "Solving the Steady State"

        sol = opt.root(self.steadystate,initss,method='lm')

        (self.k_ss,self.x_ss)=sol.x

        self.r_ss,self.w_ss,self.cy_ss,self.co_ss=self.get_ssvalues()

        if showprogress:
            print "Steady State Solved"
            print "Steady State k: ",self.k_ss
            print "Steady State x: ",self.x_ss
            print "Steady State r: ",self.r_ss
            print "Steady State w: ",self.w_ss
            print "Steady State cy: ",self.cy_ss
            print "Steady State co: ",self.co_ss

        opt_val = self.VFI(initvfi,showprogress)






        print "Now solving Value Function Iteration"



