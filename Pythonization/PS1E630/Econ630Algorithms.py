from __future__ import division
import csv
import time
import numpy as np
import scipy as sp
import Functions as FN
import scipy.optimize as opt
import scipy.interpolate as interpol
import decimal
from matplotlib import pyplot as plt
from pprint import pprint
from mpl_toolkits.mplot3d import Axes3D
from pylab import savefig


class Econ_Algorithm(object):



    def __init__(self, Parameters):

        #Loading the model parameters
        (self.A,self.alpha,self.beta,self.delta,self.eta,\
                self.N,self.tol,self.maxIter,\
                self.upper,self.lower,self.shoottol) = Parameters

        #This particular model has a closed form solution for the steady state, we
        #choose theta, so that it solves h_ss to be .33

        self.theta = (1/0.33-1)*(1-self.beta*(1-self.delta))*(1-self.alpha)/\
                (1-self.beta*(1-self.delta)-self.alpha*self.beta*self.delta)


        self.kss =  0.33/((1/self.beta-1+self.delta)/(self.alpha*self.A))**(1/(1-self.alpha))
        self.css = self.A*self.kss**self.alpha*0.33**(1-self.alpha)-self.delta*self.kss
        self.hss = .33;

    def hsolve(self,k_grid,n_pts,hsolveparams):

        hmat=np.zeros((n_pts,n_pts))
        for i in xrange(n_pts):
            hsolveparams[-1]=k_grid[i]
            h_guess=self.hss*np.ones(n_pts)

            sol = opt.root(FN.impliedhvec,h_guess,args=hsolveparams,method='lm')

            h=sol.x
            H=np.real(h)

            H[H<0]=0
            H[H>1]=1

            hmat[i,:]=H


        return hmat

    def cons(self,kgrid,hmat,i,j):
        c=self.A*kgrid[i]**self.alpha*hmat[i,j]**(1-self.alpha)-kgrid[j]\
                +(1-self.delta)*kgrid[i]

        return c

    def util(self,scale,c,hmat,i,j):
        u=scale*((c*(1-hmat[i,j])**self.theta)**(1-self.eta)-1)\
                /(1-self.eta)

        return u




    def VFI(self,progress,n_pts,guess,scale=1,reuse=None,MG=False,Mono=False,Concave=False):

        if Mono==True and Concave==True:
            print "Monotonicity and Concavity are both true, setting Concave to False"
            Concave=False

        #Evenly space grid points 
        k_grid = np.linspace(self.lower*self.kss,self.upper*self.kss,n_pts)

        initialguess=np.zeros(n_pts)

        dist = 1
        counter = 1

        v_possibilities=np.zeros((n_pts,n_pts))

        Value=np.zeros(n_pts)
        Index=np.zeros(n_pts)
        Policy=np.zeros(n_pts)


        if progress:
            print "Beginning solve for labor"

        if reuse==None:
            hsolveparams=np.zeros(n_pts+6)
            hsolveparams[0]=n_pts
            hsolveparams[1:n_pts+1]=k_grid
            hsolveparams[-2] = self.theta
            hsolveparams[-3] = self.delta
            hsolveparams[-4] = self.alpha
            hsolveparams[-5] = self.A

            hmat=self.hsolve(k_grid,n_pts,hsolveparams)

        else:
            hmat=reuse

        #print hmat
        cmat=np.zeros((n_pts,n_pts))
        umat=np.zeros((n_pts,n_pts))
        v_old=initialguess
        print "Starting VFI"

        if Concave!= True:
            while dist>self.tol:
                for i in xrange(n_pts):
                    h = hmat[i,:]
                    cmat[i,:]=self.A*k_grid[i]**self.alpha*h**(1-self.alpha)-k_grid[:]+\
                            (1-self.delta)*k_grid[i]
                    umat[i,:]=scale*((cmat[i,:]*(1-h)**self.theta)**(1-self.eta)-1)\
                            /(1-self.eta)

                    umat[i,:][umat[i,:]==np.inf]=-np.inf

                    umat[i,:][umat[i,:]==np.nan]=-np.inf

                    cmat[i,:][cmat[i,:]<0]=-np.inf
                
                    v_possibilities[i,:]=umat[i,:]+self.beta*v_old[:]
            
                if Mono!=True:
                    Value[:]=np.max(v_possibilities,axis=1)
                    Index[:]=np.argmax(v_possibilities,axis=1)

                else:
                    for i in xrange(n_pts):
                        for j in xrange(n_pts):
                            if v_possibilities[i,j]>v_possibilities[i,j+1]:
                                Value[i]=v_possibilities[i,j]
                                Index[i]=j
                                break

                            else:
                                pass
        
                dist = np.linalg.norm(Value-v_old)

                if progress:
                    print "The distance is: ",dist," with iteration number: ",counter

                if dist>self.tol:
                    v_old=Value
                    counter+=1

            if counter==self.maxIter:
                raise Exception("Took too long to converge")

        else:
            while dist>self.tol:
                t=0
                for i in xrange(n_pts):
                    j_min = t
                    j_max = n_pts-1
                    #print i
                    while j_max-j_min>2:
                        #print j_max-j_min
                        j_u=np.floor((j_min+j_max)/2)
                        j_v = j_u+1

                        c_u = self.cons(k_grid,hmat,i,j_u)

                        c_v = self.cons(k_grid,hmat,i,j_v)

                        if c_u<0:
                            v_possibilities[i,j_u]=-10**8-j_u
                        else:
                            u=self.util(scale,c_u,hmat,i,j_u)

                            v_possibilities[i,j_u]=np.real(u+self.beta*v_old[j_u])

                        if c_v<0:
                            v_possibilities[i,j_v]=-10**8-j_v
                        else:
                            u = self.util(scale,c_v,hmat,i,j_v)
                            v_possibilities[i,j_v]=np.real(u+self.beta*v_old[j_v])

                        if v_possibilities[i,j_u]<v_possibilities[i,j_v]:
                            j_min=j_u
                        else:
                            j_max=j_v

                    if j_max-j_min==2:
                        c_jmin=self.cons(k_grid,hmat,i,j_min)
                        c_jmin1=self.cons(k_grid,hmat,i,j_min+1)
                        c_jmax=self.cons(k_grid,hmat,i,j_max)

                        if c_jmin<0:
                            v_possibilities[i,j_min]=-10**8-j_min
                        else:
                            u = self.util(scale,c_jmin,hmat,i,j_min)
                            v_possibilities[i,j_min] = u+self.beta*v_old[j_min]

                        if c_jmin1<0:
                            v_possibilities[i,j_min+1]=-10**8-j_min+1
                        else:
                            u = self.util(scale,c_jmin1,hmat,i,j_min+1)
                            v_possibilities[i,j_min+1] = u+self.beta*v_old[j_min+1]

                        if c_jmax<0:
                            v_possibilities[i,j_max]=-10**8-j_max
                        else:
                            u = self.util(scale,c_jmax,hmat,i,j_max)
                            v_possibilities[i,j_max] = u+self.beta*v_old[j_max]

                        Value[i] = np.max([v_possibilities[i,j_min],v_possibilities[i,j_min+1],\
                                v_possibilities[i,j_max]])

                        t = np.argmax([v_possibilities[i,j_min],v_possibilities[i,j_min+1],\
                                v_possibilities[i,j_max]])

                        t += j_min-1

                        Index[i]=t
                    if j_max-j_min<2:
                        if c_jmin<0:
                            v_possibilities[i,j_min]=-10**8-j_min
                        else:
                            u = self.util(scale,c_jmin,hmat,i,j_min)
                            v_possibilities[i,j_min] = u+self.beta*v_old[j_min]

                        if c_jmax<0:
                            v_possibilities[i,j_max]=-10**8-j_max
                        else:
                            u = self.util(scale,c_jmax,hmat,i,j_max)
                            v_possibilities[i,j_max] = u+self.beta*v_old[j_max]

                        Value[i] = np.max([v_possibilities[i,j_min],v_possibilities[i,j_max]])

                        t = np.argmax([v_possibilities[i,j_min],v_possibilities[i,j_max]])

                        if t==0:
                            t=j_min
                        else:
                            t=j_max

                        Index[i]=t
        
                dist = np.linalg.norm(Value-v_old)

                if progress:
                    print "The distance is: ",dist," with iteration number: ",counter

                if dist>self.tol:
                    v_old=Value
                    counter+=1

                if counter==self.maxIter:
                    raise Exception("Took too long to converge")


        if MG==False:
            if scale==1:
                self.nVFIkgrid=k_grid
                self.nVFIValue=Value
                self.nVFIIndex=Index
                self.nVFIHmat=hmat

                Index=Index.astype(int)
                self.nVFIPolicy=k_grid[Index]

            else:
                self.NSSVFIkgrid=k_grid
                self.NSSVFIValue=Value
                self.NSSVFIIndex=Index

                Index=Index.astype(int)
                self.NSSVFIPolicy=k_grid[Index]

        else:
            self.MGVFIkgrid=k_grid
            self.MGVFIvalue=Value
            self.MGVFIHmat=hmat
            self.MGVFIIndex=Index
            Index=Index.astype(int)
            self.MGVFIPolicy=k_grid[Index]

    def PFI(self,initguess,n_pts,progress,scale=1,poor_n=15,reuse=None,poormans=False):

        #Evenly space grid points 
        k_grid = np.linspace(self.lower*self.kss,self.upper*self.kss,n_pts)

        dist = 1
        counter = 1

        if progress:
            print "Beginning solve for labor"

        if reuse==None:
            hsolveparams=np.zeros(n_pts+6)
            hsolveparams[0]=n_pts
            hsolveparams[1:n_pts+1]=k_grid
            hsolveparams[-2] = self.theta
            hsolveparams[-3] = self.delta
            hsolveparams[-4] = self.alpha
            hsolveparams[-5] = self.A

            hmat=self.hsolve(k_grid,n_pts,hsolveparams)

        else:
            hmat=reuse

        cmat=np.zeros((n_pts,n_pts))
        umat=np.zeros((n_pts,n_pts))

        for i in xrange(n_pts):
            h = hmat[i,:]
            cmat[i,:]=self.A*k_grid[i]**self.alpha*h**(1-self.alpha)-k_grid[:]+\
                    (1-self.delta)*k_grid[i]
            umat[i,:]=scale*((cmat[i,:]*(1-h)**self.theta)**(1-self.eta)-1)\
                            /(1-self.eta)

            umat[i,:][umat[i,:]==np.inf]=-10e9

            umat[i,:][umat[i,:]==np.nan]=-10e9

            cmat[i,:][cmat[i,:]<0]=-10e9

        dist=1
        iterate=0
        v_0=np.reshape(initguess,(1,n_pts))

        if poormans==False:
            while dist>self.tol:
                val=umat+np.dot(self.beta*np.ones((n_pts,1)),v_0)
                index = np.argmax(val.T,axis=0)

                if iterate>self.maxIter:
                    raise Exception("Took too many iterations!")

                Q = np.zeros((n_pts,n_pts))
                stepsize=n_pts*(n_pts-1)
                j = index + np.arange(n_pts,step=stepsize)
                Q[j] = 1
                Q=Q.T

                ut = umat.T
                ut=np.reshape(ut,(n_pts*n_pts))
                u=ut[j]
                u=u.T

                v_prime=np.dot(np.linalg.inv(np.identity(n_pts)-self.beta*Q),u)

                dist=np.linalg.norm(v_prime-v_0)

                v_0=np.reshape(v_prime,(1,n_pts))

                iterate+=1
                if progress:
                    print "Iteration Number: ",iterate," Distance: ",dist

        else:
            while dist>self.tol:
                val=umat+np.dot(self.beta*np.ones((n_pts,1)),v_0)
                index = np.argmax(val.T,axis=0)

                if iterate>self.maxIter:
                    raise Exception("Took too many iterations!")

                Q = np.zeros((n_pts,n_pts))
                stepsize=n_pts*(n_pts-1)
                j = index + np.arange(n_pts,step=stepsize)
                Q[j] = 1
                Q=Q.T

                ut = umat.T
                ut=np.reshape(ut,(n_pts*n_pts))
                u=ut[j]
                u=np.reshape(u,(300,1))
                
                newv=np.reshape(v_0,(n_pts,1))
                for t in xrange(poor_n):
                    newv=u+self.beta*np.dot(Q,newv)

                v_prime=np.reshape(newv,(1,n_pts))


                dist = np.linalg.norm(v_prime-v_0)

                v_0=v_prime
                iterate+=1
                if progress:
                    print "Iteration Number: ",iterate," Distance: ",dist


        return None

    def ShootingAlgorithm(self,A_trans,T,progress):
        decimal.getcontext().prec=64

        k_ss=np.zeros(2)
        c_ss=np.copy(k_ss)

        for i in xrange(2):
            k_ss[i]=((1/self.beta+self.delta-1)/(A_trans[i]*self.alpha))**(1/(self.alpha-1))
            c_ss[i]=A_trans[i]*k_ss[i]**self.alpha-self.delta*k_ss[i]

        dist=1
        c0_min=0
        c0_max=7
        iterate =0

        c0=(c0_min+c0_max)/2
        while dist>self.shoottol:

            c=np.zeros(T)
            k=np.zeros(T)

            k[0]=A_trans[-1]*k_ss[0]**self.alpha+(1-self.delta)*k_ss[0]-c0
            c[0]=c0*self.beta*(A_trans[-1]*self.alpha*k[0]**(self.alpha-1)+1-self.delta)

            for i in xrange(1,T):
                k[i] = np.real(decimal.Decimal(A_trans[-1]*k[i-1]**self.alpha+(1-self.delta)*\
                        k[i-1]-c[i-1]))
                c[i] = np.real(decimal.Decimal(c[i-1]*self.beta*(A_trans[-1]*self.alpha*\
                        k[i]**(self.alpha-1)+1-self.delta)))

                if k[i]<0 or c[i]<0:
                    break

            if k[i]<k_ss[-1]:
                c0_max=c0
            else:
                c0_min=c0

            c0=(c0_max+c0_min)/2

            dist = abs(k[i]-k_ss[-1])
            iterate+=1

            if progress:
                print "Iteration Number: ",iterate," Distance: ",dist

        self.cshootpath=c
        self.kshootpath=k





