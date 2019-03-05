# Name: Hongchao Pan
# Date: 7/18/2016
# This file contains code for numerical methods, including: Bisection method,
# Newton's methond, Secant method

from __future__ import division
# Get the reasonable approximation of x/y (true division), x//y(floor division)
import math

# Given fx
def f_x(x):
#    return math.pow(x,4)-5*x*x+4-1/(1+math.exp(math.pow(x,3)))
# For HW Q6 bootstrapping in 2016 summer refresher seminar

    r005=0.0506356160
    r01=0.0493696003
    r015=0.0475566012
    r02=0.0457436021
    r025=0.0439306021
    r03=0.0421176038
    r035=(0.5*x+1.5*r03)/2
    r04=(r03+x)/2
    r045=(1.5*x+0.5*r03)/2

    #print "r005:",r005
    #print "r01:",r01,
    #print "r015: ",r015
    #print"r02: ",r02
    #print"r025: ",r025
    #print "x:",x
    return 104-3*math.exp(-0.5*r005)-3*math.exp(-r01)-3*math.exp(-1.5*r015)- \
           3*math.exp(-2*r02)-3*math.exp(-2.5*r025)-3*math.exp(-3*r03)- \
           3*math.exp(-3.5*r035)-3*math.exp(-4*r04)-3*math.exp(-4.5*r045)-103*math.exp(-5*x)




def f_prime(x):
    #return 4*math.pow(x,3)-10*x+3*x*x*math.exp(-x*x*x)/math.pow((1+math.exp(-x*x*x)),2)
# For HW Q6 bootstrapping in 2016 summer refresher seminar
    r005=0.0506356160
    r01=0.0493696003
    r015=0.0475566012
    r02=0.0457436021
    r025=0.0439306021
    r03=0.0421176038
    r035=(0.5*x+1.5*r03)/2
    r04=(r03+x)/2
    r045=(1.5*x+0.5*r03)/2
    return 3.5*0.5*3*0.5*math.exp(-3.5*r035) + 2*3*math.exp(-4*r04) + \
           3*4.5*1.5*0.5*math.exp(-4.5*r045) + 5*103*math.exp(-5*x)


def bisection_method(a,b,f,tol_int,tol_approx):
    x_left=a
    x_right=b
    x_result=0

    while(max(abs(f(x_left)),abs(f(x_right)))>tol_approx) or ((x_right-x_left)>tol_int):
        x_result=(x_left+x_right)/2
        if (f(x_left)*f(x_result))<0:
            x_right=x_result    # active interval [x_left, x_result]
        else:
            x_left=x_result     # active interval [x_result, x_right]

    return x_result

def newtons_method(x0,f,tol_approx,tol_consec):
    x_result=x0
    x_old=x0-1

    #while (abs(f(x_result))>tol_approx) or (abs(x_result-x_old)>tol_consec):
    while(abs(x_result-x_old)>tol_consec):
        x_old=x_result
        x_result=x_old-f(x_old)/f_prime(x_old)
        #print "x_old: ",x_old
        print "Updating result: ",x_result

    return x_result

def secant_method(x00,x0,f,tol_approx,tol_consec):
    x_result=x0
    x_old=x00
    x_oldest=0

    while(abs(f(x_result))>tol_approx) or (abs(x_result-x_old)>tol_consec):
        x_oldest=x_old
        x_old=x_result
        x_result=x_old-f(x_old)*(x_old-x_oldest)/(f(x_old)-f(x_oldest))

    return x_result

if __name__ == "__main__":
    # Test the example in textbook P137-P143
    #print bisection_method(-2,3,f_x,1e-6,1e-9)
    #print newtons_method(-3,f_x,1e-9,1e-6)
    #print newtons_method(-0.5,f_x,1e-9,1e-6)
    #print newtons_method(3,f_x,1e-9,1e-6)
    #print newtons_method(0.5,f_x,1e-9,1e-6)
    #print secant_method(-3.01,-3,f_x,1e-9,1e-6)

# For HW Q6 in 2016 summer refresher seminar
    print newtons_method(0.05,f_x,1e-9,1e-6)
    #result=newtons_method(0.05,f_x,1e-9,1e-6)
   # print (result-f_x(result)/f_prime(result))

    # For HW Q6
    r00=0.05
    r005=0.0506356160
    r01=0.0493696003
    r015=0.0475566012
    r02=0.0457436021
    r025=0.0439306021
    r03=0.0421176038
    r035=0.0442946011
    r04=0.0464715983
    r045=0.0486485955
    r05=0.0508255928
    print 0.5*r05-0.5*r03










