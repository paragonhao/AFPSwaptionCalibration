#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:59:00 2019

@author: paragonhao
"""


from scipy.stats import norm
# norm.cdf(0)
import numpy as np
from swptnG2PPAF import SWPTNG2PPAF
from utils import Utils
import math
from scipy.optimize import fsolve
import datetime


# alpha, beta, sigma, eta, rho
g2params = [2.8187,0.035,0.0579,0.0091,-0.999]

# parametmers of the swaption to calibrate  
# tenor of the IRS
tenor = 1

# total notional value for the IRS
notional = 1000

# fixed rate
X = 0.02

# determine if it is a payer or receiver swaption
isPayer = 1

# the maturity of the option on the swap
maturity = 2

# frequency of the payment 
payFreq = 0.5

#load DFfactors
dffactors = Utils.readDF('DFFactor.csv')


def swaptionPrice(g2params, tenor, maturity, notional, X, isPayer):
    
    # five params in an array 
    alpha = g2params[0]
    beta = g2params[1]
    sigma = g2params[2]
    eta = g2params[3]
    rho = g2params[4]
    
    P_0_T = Utils.getZCB_0_T("07/05/19", maturity, dffactors)
    
    discountedNotional = notional * isPayer * P_0_T
    
    sigma_x = SWPTNG2PPAF.sigma_x(sigma, alpha, tenor)
    
    sigma_y = SWPTNG2PPAF.sigma_y(eta, beta, tenor)
    
    mu_x = SWPTNG2PPAF.Mu_x(g2params, tenor)
    
    mu_y = SWPTNG2PPAF.Mu_y(g2params, tenor)
    
    rho_xy = SWPTNG2PPAF.rho_xy(g2params, sigma_x, sigma_y, tenor)
    
    # first part, constant before the integral
    temp1 = notional * P_0_T
    
    # get the variable c 
    paymentTimes = tenor/payFreq
    yearFractions = np.repeat(payFreq, paymentTimes, axis=0)
    c_i = yearFractions * X
    c_i[-1] = c_i[-1] + 1
    
    t_i = pd.DataFrame(np.arange(maturity+payFreq, tenor+maturity + payFreq, payFreq))
    
    t_i_df = Utils.getDFForT_is(dffactors, tenor, maturity, payFreq)
    
    
    # numerically find the solution
    interval = 0.01
    upper_bound = mu_x + 10 * sigma_x
    lower_bound = mu_x - 10 * sigma_x
    xList = np.arange(lower_bound, upper_bound, interval)
    
    swptnPrice = 0 
    
    for x in xList:
        pdf_value = math.exp(-0.5 * (((x-mu_x)/sigma_x) ** 2)) / (sigma_x * math.sqrt( 2 * math.pi))
        
        def y_bar_solver(y_bar):
            sum_all = 0 
            for i in range(len(t_i)):
                sum_all += c_i[i] * SWPTNG2PPAF.ZCB_Func_A(g2params, t_i_df.iloc[i], P_0_T, t_i.iloc[i], tenor) * \
                math.exp(-x * SWPTNG2PPAF.B1Formula(alpha, tenor, t_i.iloc[i])) * math.exp( -y_bar * SWPTNG2PPAF.B2Formula(beta, tenor, t_i.iloc[i])) 
            return sum_all
        
        y_bar = fsolve(y_bar_solver, 1)
        print(y_bar)
        
        sqrt_rho = math.sqrt(1 - rho_xy ** 2)
        
        h_1_x = ((y_bar - mu_y)/(sigma_y * sqrt_rho)) - ((rho_xy * (x - mu_x)) / (sigma_x * sqrt_rho))
        
        temp2 = norm.cdf(-1 * isPayer * h_1_x)
        
        temp3 = 0 
        
        for j in range(len(t_i)):
            
            currB2 = SWPTNG2PPAF.B2Formula(beta, tenor, t_i.iloc[j])
            
            lambda_i = c_i[j] * SWPTNG2PPAF.ZCB_Func_A(g2params, t_i_df.iloc[j], P_0_T, t_i.iloc[j],tenor) * math.exp(-x * SWPTNG2PPAF.B1Formula(alpha, tenor, t_i.iloc[j]))
            
            kappa_i = -currB2 * (mu_y - (currB2 * 0.5 * (sigma_y**2) * (1 - rho_xy ** 2)) + (rho_xy * sigma_y * (x - mu_x)/sigma_x))
            
            h_2_x = h_1_x + currB2 * sigma_y * sqrt_rho
            
            temp3 += lambda_i * math.exp(kappa_i) * norm.cdf(-h_2_x * isPayer)
        
        swptnPrice += (temp2 - temp3) * pdf_value
        
    print("P_0_T: {}".format(P_0_T))
    print("integral Value: {}".format(swptnPrice))
    return swptnPrice * notional * P_0_T


swaptionPrice(g2params,tenor, maturity, notional, X, isPayer)















