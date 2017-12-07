import numpy as np
import matplotlib.pyplot as plt


def impliedhvec(h,params):

    N=params[0]
    kcontrol=params[1:N+1]
    (A,alpha,delta,theta,k_state)=params[N+1:]


    k=k_state
    kprime=kcontrol


    k_vec = np.ones(N)*k

    resid=theta*(A*k_vec**alpha*h**(1-alpha)+(1-delta)*k_vec-kprime)-(1-h)*(1-alpha)\
            *A*(k_vec/h)**alpha

    resid=np.real(resid)

    return resid


def plotVarPol(grid,value,policy):

    #Plotting the Value Function
    plt.title("Value Function")
    plt.plot(grid,value)
    plt.xlabel("Capital Grid")
    plt.ylabel("Value Function")
    plt.show()

    #Plotting the Policy Function
    plt.title("Policy Function")
    plt.plot(grid,policy)
    plt.xlabel("Capital Grid")
    plt.ylabel("Policy Function")
    plt.show()

def CalcResid(kgrid,value,index,hgrid,params):

    (A, alpha, beta, delta, eta,theta,n)=params

    residual = np.zeros(n)


    for i in xrange(n):
        c_rule = A*kgrid[i]**alpha*hgrid[i,index[i]]**(1-alpha)+(1-delta)*kgrid[i]-kgrid[index[i]]
        c_prime_rule = A*kgrid[index[i]]**alpha*hgrid[index[i],index[index[i]]]**(1-alpha)+(1-delta)*\
                kgrid[index[i]]-kgrid[index[index[i]]]

        rk_prime = alpha*A*kgrid[index[i]]**(alpha-1)*hgrid[index[i],index[index[i]]]**(1-alpha)

        residual[i] = 1-beta*(c_prime_rule*(1-hgrid[index[i],index[index[i]]]))**theta**(-eta)*\
                (1-hgrid[index[i],index[index[i]]])**theta*(rk_prime+1-delta)/((c_rule*(1-hgrid[i,index[i]])**theta)**(-eta)\
                *(1-hgrid[i,index[i]])**theta)


    mean_residual = np.mean(residual)
    max_residual = np.max(residual)

    print "Mean Residual: ", mean_residual
    print "Max Residual: ", max_residual



    return mean_residual, max_residual

def PlotTransPath(T,c_path,k_path):

    time = np.arange(T)

    plt.subplot(121)
    plt.plot(time,c_path)
    plt.xlabel("Time")
    plt.ylabel("Consumption Value")
    plt.title("Transition Path Consumption")

    plt.subplot(122)
    plt.plot(time,k_path)
    plt.xlabel("Time")
    plt.ylabel("Capital Value")
    plt.title("Transition Path Capital")

    plt.show()

