#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:27:24 2019

@author: paragonhao
"""

from scipy.stats import norm
# norm.cdf(0)
import numpy as np
from swptnG2PPAF import SWPTNG2PPAF
import math
from scipy.optimize import fsolve
import datetime

# alpha, beta, eta, sigma, rho
g2params = [2.8187,0.035,0.0579,0.0091,-0.999]

# parametmers of the swaption to calibrate  
tenor = 1
notional = 100 
X = 0.02
isPayer = 1
maturity = 2
payFreq = 0.5

# read in Discount factor
DFactors = pd.read_csv('DFFactor.csv')
DFactors['Date'] = pd.to_datetime(DFactors.Date)

# find the date of the maturity for the swap
startDate = "07/05/19"
startDate = datetime.datetime.strptime(startDate, "%m/%d/%y")
numDays = math.ceil(tenor * 365)

# find out the time and discount factor 
P_0_T = DFactors[DFactors['Term[day]'] < numDays].iloc[-1]['DF.mid']
t_i = pd.DataFrame(np.arange(tenor+payFreq, tenor+maturity + payFreq, payFreq))
t_i_df = t_i.apply(lambda x: getDF(x), axis= 1)

# function to get the DF function based on the number of days
def getDF(t):
    numDays = math.ceil(t * 365)
    return DFactors[DFactors['Term[day]'] < numDays].iloc[-1]['DF.mid']
    

def swaptionPrice(g2params, tenor, maturity, notional, X, isPayer):
    
    interval = 0.001
    
    # get the variable c 
    paymentTimes = maturity/payFreq
    yearFractions = np.repeat(payFreq, paymentTimes, axis=0)
    c_i = yearFractions * X
    c_i[-1] = c_i[-1] + 1
    
    # first part of the equation outside the integral;
    discountedNotional = notional * isPayer * P_0_T
    
    sigma_x = GA.sigma_x(g2params, tenor)
    sigma_y= GA.sigma_y(g2params, tenor)
    
    mu_x = -1 * GA.Moment_T_x(g2params, 0, tenor, tenor)
    mu_y = -1 * GA.Moment_T_y(g2params, 0, tenor, tenor)
    
    upper_bound = mu_x + 10 * sigma_x
    lower_bound = mu_x - 10 * sigma_x
    
    rho_xy = GA.rho_xy(g2params, tenor, sigma_x, sigma_y)
    
    # numerically solve the integration 
    xList = np.arange(lower_bound, upper_bound, interval)
    
    price = 0 
    
    for x in xList:
        
        # first part in the integral
        temp1 = math.exp(-0.5 * (((x - mu_x)/sigma_x) ** 2)) / (sigma_x * math.sqrt( 2 * math.pi))
        
        def y_bar_solver(y):
            sum_all = 0 
            for i in range(len(t_i)):
                sum_all += c_i[i] * GA.ZCB_Func_A(g2params, P_0_T, t_i_df.iloc[i], tenor, t_i.iloc[i]) * math.exp(-1 * GA.ZCB_Func_B(g2params.alpha, tenor, t_i.iloc[i]) * x - GA.ZCB_Func_B(g2params.beta, tenor, t_i.iloc[i]) * y)
                return sum_all
        
        y_bar = fsolve(y_bar_solver, 1)
        
        print(y_bar)
        
        h_1_x = ((y_bar - mu_y)/(sigma_y * math.sqrt(1 - rho_xy ** 2))) - ((rho_xy * (x - mu_x))/(sigma_x * math.sqrt(1 - rho_xy ** 2)))
        
        temp2 = norm.cdf(-1 * isPayer * h_1_x)
        
        temp3 = 0
        
        for i in range(len(t_i)):
            lambda_i = c_i[i] * GA.ZCB_Func_A(g2params, P_0_T, t_i_df.iloc[i], tenor, t_i.iloc[i]) * math.exp(-1 * GA.ZCB_Func_B(g2params.alpha, tenor, t_i.iloc[i]) * x)
            kappa_i = -1 * GA.ZCB_Func_B(g2params.beta, tenor, t_i.iloc[i]) * (mu_y - \
                                        0.5 * (1 - rho_xy ** 2) * (sigma_y ** 2) * GA.ZCB_Func_B(g2params.beta, tenor, t_i.iloc[i]) + \
                                        rho_xy * sigma_y * (x - mu_x)/sigma_x)
            
            h_2_x = h_1_x + GA.ZCB_Func_B(g2params.beta, tenor, t_i.iloc[i]) * sigma_y * math.sqrt(1 - rho_xy ** 2)
            
            temp3 += lambda_i * math.exp(kappa_i) * norm.cdf(-1 * isPayer * h_2_x)

        price += temp1 * (temp2 - temp3) * interval
    
    return price 
        

swaptionPrice(g2params, tenor, maturity, notional, X, isPayer)
        
        
swaptionPricer()
        
        
        
        
    
    