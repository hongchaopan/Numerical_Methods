# This file is written for MTH9821 homework 4, Part 4
# Hongchao Pan All right (c) reserved 2016

from __future__ import absolute_import, division, print_function
import math
import numpy as np
import black_scholes as BS
import Binomial_European as BE  # European options
import Binomial as BA           # American options
import csv

def var_reduction(t,S,K,T,sigma,q,r,N,option_type=None,method_type=None):
    '''
    Calculate the value of options via variance reduction method
    :param t:       start time
    :param S:       Spot price
    :param K:       Strike price
    :param T:       Maturity
    :param sigma:   volatility
    :param q:       dividend rate
    :param r:       risk-free rate
    :param N:       steps of binomial tree
    :param option_type: CALL or PUT
    :param method_type: Binomial, Average Binomial, BBS, BBSR
    :return:        value of options via variance reduction method
    '''
    if method_type is None:
        print("No methods selected. Please choose one of the following: Bino, Average_Bino, BBS, BBSR")
    elif method_type.upper()=="BINO":
        V,DELTA,GAMMA,THETA=BA.Binomial_American(S,K,T,sigma,q,r,N,option_type)+\
          BS.black_scholes(t,S,K,T,sigma,r,q,option_type)-\
          BE.Binomial_European(S, K, T, sigma, q, r, N, option_type)
    elif method_type.upper()=="AVERAGE_BINO":
        V, DELTA, GAMMA, THETA=BA.Average_binomial_American(S,K,T,sigma,q,r,N,option_type)+\
          BS.black_scholes(t,S,K,T,sigma,r,q,option_type)-\
          BE.Average_binomial_European(S,K,T,sigma,q,r,N,option_type)
    elif method_type.upper()=="BBS":
        V, DELTA, GAMMA, THETA=BA.BBS_American(t,S,K,T,sigma,q,r,N,option_type)+\
          BS.black_scholes(t,S,K,T,sigma,r,q,option_type)-\
          BE.BBS_European(t,S,K,T,sigma,q,r,N,option_type)
    elif method_type.upper()=="BBSR":
        V, DELTA, GAMMA, THETA=BA.BBSR_American(t,S,K,T,sigma,q,r,N,option_type)+\
          BS.black_scholes(t,S,K,T,sigma,r,q,option_type)-\
          BE.BBSR_European(t,S,K,T,sigma,q,r,N,option_type)

    return V,DELTA,GAMMA,THETA

def error(V_approx, V_exact):
    '''
    Calculate the error between approximation value and exact value
    :param V_approx:
    :param V_exact:
    :return:
    '''
    n=np.size(V_approx) # Get the size of elements
    error=np.zeros(n)
    for i in range(n):
        error[i]=abs(V_approx[i]-V_exact[i])

    return error

def error_linear(V_approx,V_exact,steps):
    n = np.size(V_approx)  # Get the size of elements
    error = np.zeros(n)
    for i in range(n):
        error[i] = steps[i]*abs(V_approx[i] - V_exact[i])

    return error

def error_qudratic(V_approx,V_exact,steps):
    n = np.size(V_approx)  # Get the size of elements
    error = np.zeros(n)
    for i in range(n):
        error[i] = steps[i] * steps[i] * abs(V_approx[i] - V_exact[i])

    return error

def get_step(n):

    step=np.zeros(n)
    step[0]=int(10)
    for i in range(1,n):
        step[i]=int(2*step[i-1])

    return step



if __name__ =="__main__":
    nn=8        # number of Ns
    steps=get_step(nn)
    print (steps)

    #Parameters
    K=40;S=41;q=1/100;sigma=30/100;r=3/100; T=1; t=0
    N=10001 # Steps of binomial tree

    # Get the exact value
    V_exact,delta_exact,gamma_exact,theta_exact=BA.Average_binomial_American(S,K,T,sigma,q,r,N,"PUT")
    print("The exact value is: ", V_exact)

    # Get the value of Binomial tree
    V_bino=np.zeros(nn)
    V_var_bino=np.zeros(nn)
    delta_bino=np.zeros(nn)
    gamma_bino=np.zeros(nn)
    theta_bino=np.zeros(nn)
    delta_var_bino=np.zeros(nn)
    gamma_var_bino=np.zeros(nn)
    theta_var_bino=np.zeros(nn)

    # Get the value of Average binomial
    V_avg_bino = np.zeros(nn)
    V_var_avg = np.zeros(nn)
    delta_avg = np.zeros(nn)
    gamma_avg = np.zeros(nn)
    theta_avg = np.zeros(nn)
    delta_var_avg = np.zeros(nn)
    gamma_var_avg = np.zeros(nn)
    theta_var_avg = np.zeros(nn)
    # Get the value of BBS
    V_BBS = np.zeros(nn)
    V_var_BBS = np.zeros(nn)
    delta_BBS = np.zeros(nn)
    gamma_BBS = np.zeros(nn)
    theta_BBS = np.zeros(nn)
    delta_var_BBS = np.zeros(nn)
    gamma_var_BBS = np.zeros(nn)
    theta_var_BBS = np.zeros(nn)
    # Get the value of BBSR
    V_BBSR=np.zeros(nn)
    V_var_BBSR=np.zeros(nn)
    delta_BBSR = np.zeros(nn)
    gamma_BBSR = np.zeros(nn)
    theta_BBSR = np.zeros(nn)
    delta_var_BBSR = np.zeros(nn)
    gamma_var_BBSR = np.zeros(nn)
    theta_var_BBSR = np.zeros(nn)

    # Get the V_exact_vec
    V_exact_vec=np.zeros(nn)
    delta_exact_vec=np.zeros(nn)
    gamma_exact_vec = np.zeros(nn)
    theta_exact_vec = np.zeros(nn)

    for i in range(nn):
        V_bino[i],delta_bino[i],gamma_bino[i],theta_bino[i]=BA.Binomial_American(S,K,T,sigma,q,r,steps[i],"PUT")
        V_var_bino[i],delta_var_bino[i],gamma_var_bino[i],theta_var_bino[i]=var_reduction(t,S,K,T,sigma,q,r,steps[i],"PUT","BINO")

        V_avg_bino[i],delta_avg[i],gamma_avg[i],theta_avg[i] = BA.Average_binomial_American(S, K, T, sigma, q, r, steps[i], "PUT")
        V_var_avg[i],delta_var_avg[i],gamma_var_avg[i],theta_var_avg[i] = var_reduction(S, K, T, sigma, q, r, steps[i], "PUT", "AVERAGE_BINO")

        V_BBS[i],delta_BBS[i],gamma_BBS[i],theta_BBS[i]=BA.BBS_American(t,S,K,T,sigma,q,r,steps[i],"PUT")
        V_var_BBS[i]=var_reduction(t,S,K,T,sigma,q,r,steps[i],"PUT","BBS")

        V_BBSR[i]=BA.BBSR_American(t,S,K,T,sigma,q,r,steps[i],"PUT")
        V_var_BBSR[i],delta_var_BBS[i],gamma_var_BBS[i],theta_var_BBS[i]=var_reduction(t,S,K,T,sigma,q,r,steps[i],"PUT","BBSR")

        V_exact_vec[i]=V_exact
        delta_exact_vec[i]=delta_exact
        gamma_exact_vec[i]=gamma_exact
        theta_exact_vec[i]=theta_exact

    # Get the errors

    # Get the value error
    error_bino=error(V_bino,V_exact_vec)
    error_bino_linear=error_linear(V_bino,V_exact_vec,steps)
    error_bino_quadratic=error_qudratic(V_bino,V_exact_vec,steps)

    error_delta_bino=error(delta_bino,delta_exact)
    error_delta_bino_linear=error_linear(delta_bino,delta_exact,steps)
    error_delta_bino_quadratic=error_qudratic(delta_bino,delta_exact,steps)

    error_gamma_bino = error(gamma_bino, gamma_exact)
    error_gamma_bino_linear = error_linear(gamma_bino, gamma_exact, steps)
    error_gamma_bino_quadratic = error_qudratic(gamma_bino, gamma_exact, steps)

    error_theta_bino=error(theta_bino,theta_exact)
    error_theta_bino_linear = error_linear(theta_bino, theta_exact,steps)
    error_theta_bino_quadratic = error_qudratic(theta_bino, theta_exact,steps)

    error_avg=error(V_avg_bino,V_exact_vec)
    error_avg_linear=error_linear(V_avg_bino,V_exact_vec,steps)
    error_avg_quadratic=error_qudratic(V_avg_bino,V_exact_vec,steps)

    error_delta_avg = error(delta_avg, delta_exact)
    error_delta_avg_linear = error_linear(delta_avg, delta_exact, steps)
    error_delta_avg_quadratic = error_qudratic(delta_avg, delta_exact, steps)

    error_gamma_avg = error(gamma_avg, gamma_exact)
    error_gamma_avg_linear = error_linear(gamma_avg, gamma_exact, steps)
    error_gamma_avg_quadratic = error_qudratic(gamma_avg, gamma_exact, steps)

    error_theta_avg = error(theta_avg, theta_exact)
    error_theta_avg_linear = error_linear(theta_avg, theta_exact, steps)
    error_theta_avg_quadratic = error_qudratic(theta_avg, theta_exact, steps)


    error_BBS=error(V_BBS,V_exact_vec)
    error_BBS_linear=error_linear(V_BBS,V_exact_vec,steps)
    error_BBS_quadratic=error_qudratic(V_BBS,V_exact_vec,steps)

    error_delta_BBS = error(delta_BBS, delta_exact)
    error_delta_BBS_linear = error_linear(delta_BBS, delta_exact, steps)
    error_delta_BBS_quadratic = error_qudratic(delta_BBS, delta_exact, steps)

    error_gamma_BBS = error(gamma_BBS, gamma_exact)
    error_gamma_BBS_linear = error_linear(gamma_BBS, gamma_exact, steps)
    error_gamma_BBS_quadratic = error_qudratic(gamma_BBS, gamma_exact, steps)

    error_theta_BBS = error(theta_BBS, theta_exact)
    error_theta_BBS_linear = error_linear(theta_BBS, theta_exact, steps)
    error_theta_BBS_quadratic = error_qudratic(theta_BBS, theta_exact, steps)

    error_BBSR=error(V_BBSR,V_exact_vec)
    error_BBSR_linear=error_linear(V_BBSR,V_exact_vec,steps)
    error_BBSR_quadratic=error_qudratic(V_BBSR,V_exact_vec,steps)

    error_delta_BBSR = error(delta_BBSR, delta_exact)
    error_delta_BBSR_linear = error_linear(delta_BBSR, delta_exact, steps)
    error_delta_BBSR_quadratic = error_qudratic(delta_BBSR, delta_exact, steps)

    error_gamma_BBSR = error(gamma_BBSR, gamma_exact)
    error_gamma_BBSR_linear = error_linear(gamma_BBSR, gamma_exact, steps)
    error_gamma_BBSR_quadratic = error_qudratic(gamma_BBSR, gamma_exact, steps)

    error_theta_BBSR = error(theta_BBSR, theta_exact)
    error_theta_BBSR_linear = error_linear(theta_BBSR, theta_exact, steps)
    error_theta_BBSR_quadratic = error_qudratic(theta_BBSR, theta_exact, steps)

    # Get the value error of variance reduction
    error_var_bino = error(V_var_bino, V_exact_vec)
    error_var_bino_linear = error_linear(V_var_bino, V_exact_vec, steps)
    error_var_bino_quadratic = error_qudratic(V_var_bino, V_exact_vec, steps)

    error_var_delta_bino = error(delta_var_bino, delta_exact)
    error_var_delta_bino_linear = error_linear(delta_var_bino, delta_exact, steps)
    error_var_delta_bino_quadratic = error_qudratic(delta_var_bino, delta_exact, steps)

    error_var_gamma_bino = error(gamma_var_bino, gamma_exact)
    error_var_gamma_bino_linear = error_linear(gamma_var_bino, gamma_exact, steps)
    error_var_gamma_bino_quadratic = error_qudratic(gamma_var_bino, gamma_exact, steps)

    error_var_theta_bino = error(theta_var_bino, theta_exact)
    error_var_theta_bino_linear = error_linear(theta_var_bino, theta_exact, steps)
    error_var_theta_bino_quadratic = error_qudratic(theta_var_bino, theta_exact, steps)

    error_var_avg = error(V_var_avg, V_exact_vec)
    error_var_avg_linear = error_linear(V_var_avg, V_exact_vec, steps)
    error_var_avg_quadratic = error_qudratic(V_var_avg, V_exact_vec, steps)

    error_var_delta_avg = error(delta_var_avg, delta_exact)
    error_var_delta_avg_linear = error_linear(delta_var_avg, delta_exact, steps)
    error_var_delta_avg_quadratic = error_qudratic(delta_var_avg, delta_exact, steps)

    error_var_gamma_avg = error(gamma_var_avg, gamma_exact)
    error_var_gamma_avg_linear = error_linear(gamma_var_avg, gamma_exact, steps)
    error_var_gamma_avg_quadratic = error_qudratic(gamma_var_avg, gamma_exact, steps)

    error_var_theta_avg = error(theta_var_avg, theta_exact)
    error_var_theta_avg_linear = error_linear(theta_var_avg, theta_exact, steps)
    error_var_theta_avg_quadratic = error_qudratic(theta_var_avg, theta_exact, steps)

    error_var_BBS = error(V_var_BBS, V_exact_vec)
    error_var_BBS_linear = error_linear(V_var_BBS, V_exact_vec, steps)
    error_var_BBS_quadratic = error_qudratic(V_var_BBS, V_exact_vec, steps)

    error_var_BBSR = error(V_var_BBSR, V_exact_vec)
    error_var_BBSR_linear = error_linear(V_var_BBSR, V_exact_vec, steps)
    error_var_BBSR_quadratic = error_qudratic(V_var_BBSR, V_exact_vec, steps)

    # Write the results to csv files
    #csvfile = file("error.csv", "wb")
    #writer = csv.writer(csvfile)
    #writer.writerows(V_bino)
    #csvfile.close()

    print ("Binomial")
    print(V_bino)
    print("******")
    print(error_bino)
    print("*******")
    print(error_bino_linear)
    print("*******")
    print(error_bino_quadratic)
    print("*******")
    print(delta_bino)
    print("*******")
    print(gamma_bino)
    print("*******")
    print(theta_bino)

    print("Average Binomial")
    print(V_avg_bino)
    print("******")
    print(error_avg)
    print("*******")
    print(error_avg_linear)
    print("*******")
    print(error_avg_quadratic)
    print("*******")
    print(delta_avg)
    print("*******")
    print(gamma_avg)
    print("*******")
    print(theta_avg)

    print("BBS")
    print(V_BBS)
    print("******")
    print(error_BBS)
    print("*******")
    print(error_BBS_linear)
    print("*******")
    print(error_BBS_quadratic)
    print("*******")
    print(delta_BBS)
    print("*******")
    print(gamma_BBS)
    print("*******")
    print(theta_BBS)

    print("BBSR")
    print(V_BBSR)
    print("******")
    print(error_BBSR)
    print("*******")
    print(error_BBSR_linear)
    print("*******")
    print(error_BBSR_quadratic)
    print("*******")
    print(delta_BBSR)
    print("*******")
    print(gamma_BBSR)
    print("*******")
    print(theta_BBSR)

    print("Variance Reduction")
    print("Binomial")
    print(V_var_bino)
    print("******")
    print(error_var_bino)
    print("*******")
    print(error_var_bino_linear)
    print("*******")
    print(error_var_bino_quadratic)
    print("*******")
    print(delta_bino)
    print("*******")
    print(gamma_bino)
    print("*******")
    print(theta_bino)

    print("Average Binomial")
    print(V_avg_bino)
    print("******")
    print(error_avg)
    print("*******")
    print(error_avg_linear)
    print("*******")
    print(error_avg_quadratic)
    print("*******")
    print(delta_avg)
    print("*******")
    print(gamma_avg)
    print("*******")
    print(theta_avg)

    print("BBS")
    print(V_BBS)
    print("******")
    print(error_BBS)
    print("*******")
    print(error_BBS_linear)
    print("*******")
    print(error_BBS_quadratic)
    print("*******")
    print(delta_BBS)
    print("*******")
    print(gamma_BBS)
    print("*******")
    print(theta_BBS)

    print("BBSR")
    print(V_BBSR)
    print("******")
    print(error_BBSR)
    print("*******")
    print(error_BBSR_linear)
    print("*******")
    print(error_BBSR_quadratic)
    print("*******")
    print(delta_BBSR)
    print("*******")
    print(gamma_BBSR)
    print("*******")
    print(theta_BBSR)


