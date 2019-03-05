# Name: Hongchao Pan
# Date: 7/17/2016
# This file contains mid-point rule, Traprzoidal rule, and Simposon's rule
# and calculate approximate value of an integral with given tolerance

from __future__ import division
# Get the reasonable approximation of x/y (true division), x//y(floor division)
import math

# Define input evaluating function
def f_x(x):
#    return math.exp(-x*x)
# For Quiz3 2016 summer refresher seminar
# For computing N(t), compute integration of e^(-0.5x^2)/sqrt(2pi)
    #return (math.exp(-0.5*x*x))/math.sqrt(2*math.pi)
    '''
    # For Q12 in refresher seminar
    t=0
    T=0.25
    S=50
    K=45
    sigma=0.25
    q=0.01
    r=0.03
    d2=(math.log(S/K)+(r-q-0.5*math.pow(sigma,2))*(T-t))/(sigma*math.sqrt(T-t))
    return math.sqrt(-math.log(x))*math.exp(-0.5*(math.log(x)+d2)*(math.log(x)+d2))/x
    '''
     # For Q13(i) in refresher seminar
    t=0
    T=3/12
    S=50
    K=50
    sigma=0.3
    q=0.02
    r=0.04
    d2=(math.log(S/K)+(r-q-0.5*math.pow(sigma,2))*(T-t))/(sigma*math.sqrt(T-t))
    # For Q13(i)
    return math.sqrt(1-x)*math.exp(-0.5*(math.log(x)/(sigma*math.sqrt(T-t))-d2)*(math.log(x)/(sigma*math.sqrt(T-t))-d2))/x
    #return math.sqrt(1/x-1)*math.exp(-0.5*(math.log(x)/(sigma*math.sqrt(T-t))+d2)*(math.log(x)/(sigma*math.sqrt(T-t))+d2))/x

# a: left endpoint; b: right endpoint; n: number of partition intervals
def mid_point_rule(a,b,f,n):
    h=(b-a)/n
    result=0
    for i in xrange(1,n+1):
        result+=f(a+(i-0.5)*h)
       # print result

    return result*h

# a: left endpoint; b: right endpoint; n: number of partition intervals
def Trapezoidal_rule(a,b,f,n):
    h= (b-a)/n
    result=(f(a)+f(b))/2

    for i in xrange(1, n):
        result += f(a + i*h)

    return result*h

# a: left endpoint; b: right endpoint; n: number of partition intervals
def simpson_rule(a,b,f,n):
    h= (b-a)/n
    result=(f(a) + f(b))/6
    #print "initial result: ",result

    for i in xrange(1, n):
        result += f(a + i*h)/3
    for i in xrange(1, n+1):
        result += 2*f(a + (i-0.5)*h)/3
        #print "Updating result: ",result*h
    return result*h

# f: mid_point/trape/simpson method, n: intervals
def approx_val_tol(f_med,a,b,f,n,tol):
    #tol=5e-7
    #print "Approximation: "
    #print "Initial n is: ", n
    result_old=f_med(a,b,f,n)
    #print "Initial result is: ", result_old
    #print "Initial result + 0.5 is: ", result_old+0.5
    n=2*n
    #print "Updating n is: ", n
    result_new=f_med(a,b,f,n)
    #print "Updating result is; ", result_new
    #print "Updating result + 0.5 is: ", result_new+0.5

    while(abs(result_old-result_new)>tol):
        result_old=result_new
        n=2*n
        result_new=f_med(a,b,f,n)
        #print "Updating n is: ", n
        #print "Updating result is: ", result_new    # To see the result in each interval
        #print "Updating result + 0.5 is: ", result_new+0.5
    return result_new


# execute program from here
if __name__=="__main__":
    #print mid_point_rule(0,2,f_x,4) # check answer with P49
    #print mid_point_rule(0,2,f_x,512)   # check answer with P49
    #print simpson_rule(0,2,f_x,4) # check answer with P49

    #print approx_val_tol(mid_point_rule,0,2,f_x,4,5e-7)
    #print approx_val_tol(simpson_rule,0,2,f_x,4,5e-7)

    #print "Below is for Quiz3 Q1 in 2016 summer refersher seminar"
    #print "Result for N(0.1): "
    #print approx_val_tol(simpson_rule,0,0.1,f_x,4,1e-12)
    #print "Result for N(0.5): "
    #print approx_val_tol(simpson_rule,0,0.5,f_x,4,1e-12)
    #print "Result for N(1.0): "
    #print approx_val_tol(simpson_rule,0,1.0,f_x,4,1e-12)
    '''
    # Q12 in refresher
    t=0
    T=0.25
    S=50
    K=45
    sigma=0.25
    q=0.01
    r=0.03
    #print approx_val_tol(simpson_rule,1e-200,1,f_x,4,1e-6)
    print approx_val_tol(simpson_rule,1e-200,1,f_x,4,1e-6)*math.exp(-r*T)*math.sqrt(sigma*math.sqrt(T)/(2*math.pi))
    #print simpson_rule(1e-300,1,f_x,4)
    #print simpson_rule(1e-300,1,f_x,4)*math.exp(-r*T)*math.sqrt(sigma*math.sqrt(T)/(2*math.pi))
    '''

    # Q13(i) in refresher
    t=0
    T=3/12
    S=50
    K=50
    sigma=0.3
    q=0.02
    r=0.04
    # For Q13 (i) and (ii)
    print approx_val_tol(simpson_rule,1e-200,1,f_x,4,1e-6)*math.exp(-r*T)*math.sqrt(K/(T*2*math.pi))/sigma
