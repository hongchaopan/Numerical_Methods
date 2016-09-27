# Name: Hongchao Pan
# Date: 7/18/2016
# This file contains code for computing black_scholes formula and its greeks
# Must include cumulative_distribution.py file

from __future__ import division, absolute_import, print_function
# Get the reasonable approximation of x/y (true division), x//y(floor division)
import math
#from cumulative_distribution import cum_dist_normal
from scipy.stats import norm
#import logging
#logging.basicConfig(format='%(module)s - %(funcName)s - %(message)s', level=logging.DEBUG)
#log = logging.getLogger(__name__)

def d_1(t,S,K,T,sigma,r,q):
    return (math.log(S/K)+(r-q+0.5*math.pow(sigma,2))*(T-t))/(sigma*math.sqrt(T-t))

def black_scholes(t,S,K,T,sigma,r,q,option_type=None):
    d1=d_1(t,S,K,T,sigma,r,q)
    d2=d1-sigma*math.sqrt(T-t)

    call_price=S*math.exp(-q*(T-t))*norm.cdf(d1)-K*math.exp(-r*(T-t))*norm.cdf(d2)
    put_price=K*math.exp(-r*(T-t))*norm.cdf(-d2)-S*math.exp(-q*(T-t))*norm.cdf(-d1)

    if option_type is None:
        print ("No option type indicated, assuming CALL.")
        #log.warning("No option type indicated, assuming CALL.")
        return call_price
    elif option_type.upper()=="CALL":
        return call_price
    elif option_type.upper()=="PUT":
        return put_price

# Define the greeks of black_scholes models
def delta_BS_call(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    return math.exp(-q*(T-t))*norm.cdf(d1)
def delta_BS_put(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    return -math.exp(-q*(T-t))*norm.cdf(-d1)

def gamma_BS(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    return math.exp(-q*(T-t))*math.exp(-0.5*d1*d1)/(S*sigma*math.sqrt(2*math.pi*(T-t)))

def vega_BS(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    return S*math.exp(-q*(T-t))*math.sqrt(T-t)*math.exp(-0.5*d1*d1)/math.sqrt(2*math.pi)

def theta_BS_call(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    d2=d1-sigma*math.sqrt(T-t)
    return -0.5*S*sigma*math.exp(-q*(T-t))*math.exp(-0.5*d1*d1)/math.sqrt(2*math.pi*(T-t))+ \
q*S*math.exp(-q*(T-t))*norm.cdf(d1)-r*K*math.exp(-r*(T-t))*norm.cdf(d2)

def theta_BS_put(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    d2=d1-sigma*math.sqrt(T-t)
    return -0.5*S*sigma*math.exp(-q*(T-t))*math.exp(-0.5*d1*d1)/math.sqrt(2*math.pi*(T-t))- \
q*S*math.exp(-q*(T-t))*norm.cdf(-d1) + r*K*math.exp(-r*(T-t))*norm.cdf(-d2)

def rho_BS_call(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    d2=d1-sigma*math.sqrt(T-t)
    return K*(T-t)*math.exp(-r*(T-t))*norm.cdf(d2)

def rho_BS_put(t,S,K,T,sigma,r,q):
    d1=d_1(t,S,K,T,sigma,r,q)
    d2=d1-sigma*math.sqrt(T-t)
    return -K*(T-t)*math.exp(-r*(T-t))*norm.cdf(-d2)

if __name__ == "__main__":
# Test call and put option
    print (black_scholes(0,42,40,0.5,0.3,0.05,0.03))
    print (black_scholes(0,42,40,0.5,0.3,0.05,0.03,"CALL"))
    print (black_scholes(0,42,40,0.5,0.3,0.05,0.03,"PUT"))
    print(delta_BS_put(0,42,40,0.5,0.3,0.05,0.03))