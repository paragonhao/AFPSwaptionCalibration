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
import math
import scipy.optimize as optimize
from normalModel import NormalModel
import datetime
import pandas as pd

# parametmers of the swaption to calibrate  
# tenor of the IRS
tenor = 1.0

# the maturity of the option on the swap
maturity = 2.0

# total notional value for the IRS
notional = 1.0

# fixed rate
fixedRate = 0.0166

# frequency of the payment 
payFreq = 0.5

#load DFfactors
termStructure = Utils.readDF('termStructure.csv')

G2 = SWPTNG2PPAF(termStructure)

def swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate):
    print("running the swaptionPricingFunction")
    isPayer = 1
    isReceiver = -1 
    
    # five params in an array 
    alpha = g2params[0]
    beta = g2params[1]
    sigma = g2params[2]
    eta = g2params[3]
    rho = g2params[4]
    print("parameters +++++++++++++++++++++++++++")
    print(g2params)
    print("+++++++++++++++++++++")
    P_0_T = G2.getTermStructure(maturity)
    
    # find out the cash flow and the payment time in terms of years
    paymentTimes = tenor/payFreq
    tau = np.repeat(payFreq, paymentTimes, axis=0)
    c_i = tau * fixedRate
    c_i[-1] = c_i[-1] + 1
    t_i = np.arange(maturity+payFreq, tenor+maturity + payFreq, payFreq)
          
    sigma_x = G2.sigma_x(sigma, alpha, maturity)
    
    sigma_y = G2.sigma_y(eta, beta, maturity)
    
    mu_x = G2.Mu_x(g2params, maturity)
    
    mu_y = G2.Mu_y(g2params, maturity)
    
    rho_xy = G2.rho_xy(g2params, sigma_x, sigma_y, maturity)
    
    txy = math.sqrt(1 - rho_xy * rho_xy)
       
    A_array = []
    Ba_array = []
    Bb_array = []
    lambda_array = np.repeat(0, len(t_i), axis=0)
    
    # define the boundary 
    interval = 0.01
    lower = mu_x - 10 * sigma_x
    upper = mu_x + 10 * sigma_x
    xList = np.arange(lower, upper, interval)
    
    for i in range(len(t_i)):

        A_array.append(G2.A(g2params, maturity, t_i[i]))

        Ba_array.append(G2.B(alpha, t_i[i] - maturity))

        Bb_array.append(G2.B(beta, t_i[i] - maturity))
    
    
    integral_result_Payer = 0 
    integral_result_Receiver = 0 
    cdf_val_1 = 0               # payer
    cdf_val_2 = 0               # receiver
    
    # numerically solves the integral
    for x in xList:
        lambda_array = []
        for i in range(len(t_i)):
            lambda_array.append(c_i[i] * A_array[i] * math.exp(-Ba_array[i] * x))

        def y_bar_solver(y):
            sum_all = 0.0
            
            for i in range(len(t_i)):
                sum_all += lambda_array[i] * math.exp(-Bb_array[i] * y)
            return sum_all - 1.0
        
        y_bar = optimize.fsolve(y_bar_solver, x0=1, xtol=1e-6)
        h1 = (y_bar - mu_y)/(sigma_y * txy) - rho_xy * (x - mu_x)/(sigma_x * txy)
        
        cdf_val_1 = norm.cdf(-isPayer * h1)
        cdf_val_2 = norm.cdf(-isReceiver * h1)
        
        try:
            for i in range(len(t_i)):

                h2 = h1 + Bb_array[i] * sigma_y * txy
                kappa = - Bb_array[i] * (mu_y - 0.5 * txy * txy * sigma_y * sigma_y * Bb_array[i] \
                                  + rho_xy * sigma_y * (x - mu_x)/sigma_x)
                cdf_val_1 -= lambda_array[i] * math.exp(kappa) * norm.cdf(-h2 * isPayer)[0]
                cdf_val_2 -= lambda_array[i] * math.exp(kappa) * norm.cdf(-h2 * isReceiver)[0]
        except OverflowError:
            print("Overflow err: ")
            
        temp_pdf = (x - mu_x) / sigma_x
        
        integral_result_Payer += math.exp(-0.5 * temp_pdf * temp_pdf) * cdf_val_1/(sigma_x * math.sqrt(2.0 * math.pi))
        integral_result_Receiver += math.exp(-0.5 * temp_pdf * temp_pdf) * cdf_val_2/(sigma_x * math.sqrt(2.0 * math.pi))
    
    # swaption price calculated using G2++ analytical solution
    swaptionPricePayer = notional * isPayer * P_0_T * integral_result_Payer
    swaptionPriceReceiver = notional * isReceiver * P_0_T * integral_result_Receiver

    # get the discount factor for each payment of the interest rate swap    
    dfs = [G2.getTermStructure(t) for t in t_i]
    
    # 0.00709 is the IV of the sum of the payer and receiver ATM swaption, so the IV for each swaption is 0.00709/2
    normModel = NormalModel(fRate = fixedRate, strike = fixedRate, impliedVol = 0.00709/2, maturity = maturity, discountFactors=dfs)


    return (swaptionPricePayer + swaptionPriceReceiver) - (normModel.getPayerSwaptionPrice() + normModel.getReceiverSwaptionPrice())
    
#g2params = [2.8187,0.035,0.0579,0.0091,-0.999]
g2params = [1,1,1,1,-0.999]
swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate)
#bnds = ((0.001, 5),(0.001, 0.05),(0.001, 0.05),(0.001, 0.05),(-0.999, 0.999))
bnds = ((0.001, 5),(0.001, 5),(0.001, 5),(0.001, 5),(-0.999, 0.999))
result = optimize.minimize(swaptionPricingFunction, g2params, args=(tenor, maturity, notional, fixedRate), bounds=bnds, method='L-BFGS-B')

if result.success:
    fitted_params = result.x
    print("current fitted value")
    print(fitted_params)
else:
    raise ValueError(result.message)
    






