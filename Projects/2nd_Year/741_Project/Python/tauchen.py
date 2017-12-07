from __future__ import division
import numpy as np
from scipy.stats import norm


def tauchen(rho,sig_e,n,disp):
    """
    Description: 
        - A python translation of a MATLAB code that finds the transition probabilities
          for a Markov chain that approximates an AR(1) in the form:

          y[t] = rho*y[t-1]+e[t], e[t] ~ N[0,sig_e^2]

          References:
          Tauchen, G. "Finite State Markov-chain approximations to
          Univariate and Vector Autogregressions," Economics Letters
          20, pp. 171-181, 1986.

          Jesus Fernandez-Villaverde
          Minneapolis 5-29-2001

          Python Translation by:
          Jeff Clawson, 11-22-2017
          University of Maryland, College Park

          Accuracy has been as closely verified as possible. Any mistakes here
          are my own.
            
    Inputs:
        - rho       = Parameter: AR(1) slope parameter
        - sig_e     = Parameter: AR(1) variance parameter
        - n         = Parameter: Number of points for the grid
        - disp      = Parameter: Amount of dispersion
    
    Other Functions Called:
        - norm.cdf      = From scipy, gives the CDF of a standard normal distribution
        - np.sqrt       = From numpy, square root
        - np.ones       = From numpy, vector of ones of a givens size
        - np.dot        = Matrix multiplication, from numpy
        - np.linalg.norm= Calculates the norm from numpy
        - np.linspace   = From numpy, calculate and evenly spaced grid

    Outputs:
        - y             = Array [n], points of approximation
        - p             = Array [n,n], transition kernel
        - pi            = Array [n], invariant distribution
        

    """


    #1. Points for y
    #End point for the grid
    stop = disp*sig_e/np.sqrt(1-rho**2)

    #By definition, the start point of the grid is the 
    #negative of the end point
    start = -stop

    #Step size for the grid, not particularly necessary for setting up the
    #grid, but it is used later
    w = 2*stop/(n-1)

    #Sets up the y grid
    y = np.linspace(start,stop,n)

    #2. Compute the Transition
    #Initiates the transition
    p = np.zeros((n,n))

    #Fills in the first column
    p_in1 = (y[0] - rho*y[:]+(w/2))/sig_e
    p[:,0] = norm.cdf(p_in1)

    #Fills in the last column
    p_inn =  (y[-1] - rho*y[:]-(w/2))/sig_e
    p[:,-1] = 1-norm.cdf(p_inn)

    #Then, cycles through and fills in the remaining columns,
    #one row at a time
    for j in xrange(n):
        #Initializes value to be plugged in for the first part
        p_injl = (y[1:-1]-rho*y[j]+(w/2))/sig_e 

        #Initializes value to be plugged in for the second part
        p_injr = (y[1:-1]-rho*y[j]-(w/2))/sig_e 

        #Fills in the row
        p[j,1:-1] = norm.cdf(p_injl)-norm.cdf(p_injr)


    #3. Computing the invariant distribution

    #Starts with a uniform distribution
    pi = np.ones(n)/n

    #Adds one to it
    piprov = pi+1

    #Initializes the distance measure
    dist = 1
    
    #Goes through and calculates a fixed point of the distribution
    while dist>.00001:
        piprov = pi
        pi = np.dot(pi,p)
        
        #Takes the norm of two vectors
        dist = np.linalg.norm(pi-piprov)

    return y, p, pi
