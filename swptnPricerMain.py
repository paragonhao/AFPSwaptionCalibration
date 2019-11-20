#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:59:00 2019

@author: paragonhao
"""
from scipy.stats import norm
import numpy as np
from swptnG2PPAF import SWPTNG2PPAF
import math
import scipy.optimize as optimize
from normalModel import NormalModel
import pandas as pd
from utils import Utils

# total notional value for the IRS
notional = 1.0

# frequency of the payment 
payFreq = 0.25

################################## data initilaisation ############################################
termStructure = Utils.readDF('data/term_structure_interploated.csv')

mktVol = Utils.readMktVolSurface('data/marketVol.csv')

optimisationGrid = Utils.readOptimGrid('data/optimisationGrid.csv')

fssGrid = Utils.readFSS('data/fssGrid.csv')

G2 = SWPTNG2PPAF(termStructure)

optimisationGrid['Vol'] = 0.0

for idx, row in optimisationGrid.iterrows():
    curVol = mktVol[row['Tenor']][row['Maturity']]
    optimisationGrid['Vol'][idx] = curVol

optimisationGrid['Fss'] = 0.0

for idx, row in optimisationGrid.iterrows():
    optimisationGrid['Fss'][idx] =  fssGrid[row['Tenor']][row['Maturity']]

######################################################################################################



################################## Optimisation Function  ############################################
def swaptionPricingOptim(g2params, optimisationGrid):
    
    rmse = 0.0
    
    for idx, row in optimisationGrid.iterrows():
        
        maturity = Utils.monthToYear(row['Maturity'])
        tenor = Utils.monthToYear(row['Tenor'])
        
        # vol is in basis points
        mktVol = row['Vol']/10000.0
        
        print("parameters +++++++++++++++++++++++++++")
        print(g2params)
        print("+++++++++++++++++++++")
        
        # fixed rate is in percentage
        fixedRate = row['Fss']/100.0
        print("Maturity: {}, tenor: {}, vol: {}, fixedRate: {}".format(maturity, tenor, mktVol, fixedRate))
        
        normP, normR, G2P, G2R = swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate, mktVol)
        currErr = (normP + normR - G2P - G2R) ** 2
        rmse += currErr
        
    print("Total RMSE: {}".format(rmse))
    return rmse/len(optimisationGrid)
######################################################################################################
     


################################## Swaption Pricing Function  ############################################
def swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate, mktVol):
#    print("running the swaptionPricingFunction")
    isPayer = 1
    isReceiver = -1 
    
    # five params in an array 
    alpha = g2params[0]
    beta = g2params[1]
    sigma = g2params[2]
    eta = g2params[3]
    rho = g2params[4]
    
    P_0_T = G2.getTermStructure(maturity)
    
    # find out the cash flow and the payment time in terms of years
    # tau: payment Time in terms of year
    # c_i: cash flow at each time of the IRS
    # t_i: time at which the payments are made
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
    interval = 0.0001
    lower = mu_x - 10 * sigma_x
    upper = mu_x + 10 * sigma_x

    xList = np.arange(lower, upper, interval)
    
    for i in range(len(t_i)): 
        A_array.append(G2.A(g2params, maturity, t_i[i]))

        Ba_array.append(G2.B(alpha, maturity, t_i[i]))

        Bb_array.append(G2.B(beta, maturity, t_i[i]))
    
    
    integral_result_Payer = 0 
    integral_result_Receiver = 0 
    cdf_val_1 = 0               # payer
    cdf_val_2 = 0               # receive
    
    # numerically solves the integral
    for x in xList:
#        print("XList max: {}, min: {}".format(max(xList), min(xList)))
        lambda_array = []
        
        for i in range(len(t_i)):
#            print("gett array values c_i: {}, A_array: {}, Ba_array: {}, x: {}".format(c_i[i], A_array[i], Ba_array[i], x))
            try:
                lambda_array.append(c_i[i] * A_array[i] * math.exp(-Ba_array[i] * x))
            except OverflowError:
                print("exponential value on Ba_array is too big {}")     
        
        # confirm if this is correct 
        def y_bar_solver(y):
            sum_all = 0.0
            
            for i in range(len(t_i)):
                sum_all += lambda_array[i] * math.exp(-Bb_array[i] * y)
            return sum_all - 1.0
        
        y_bar = optimize.fsolve(y_bar_solver, x0=1, xtol=1e-6)
        
        h1 = (y_bar - mu_y)/(sigma_y * txy) - (rho_xy * (x - mu_x))/(sigma_x * txy)
        
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
    
    normModel = NormalModel(fRate = fixedRate, strike = fixedRate, impliedVol = mktVol/2.0, maturity = maturity, discountFactors=dfs)
    
    return normModel.getPayerSwaptionPrice(), normModel.getReceiverSwaptionPrice(), swaptionPricePayer, swaptionPriceReceiver
############################################################################################################################################################



########################################### execution ########################################################
g2params = [2.8187,0.035,0.0579,0.0091,-0.999]
g2params = [5.53449534, 3.70068385, 0.001, 0.001, 0.925757890]

swaptionPricingOptim(g2params, optimisationGrid)
bnds = ((0.001, 100),(0.001, 100),(0.001, 100),(0.001, 100),(-0.999, 0.999))
result = optimize.minimize(swaptionPricingOptim, g2params, args=(optimisationGrid), bounds=bnds, method='L-BFGS-B')


if result.success:
    fitted_params = result.x
    print("current fitted value")
    print(fitted_params)
else:
    raise ValueError(result.message)
########################################################################################
 
    

########################################### compare results ########################################################
#def swaptionModelVsMkt(g2params, optimisationGrid):
#    
#    rmse = 0.0
#    
#    for idx, row in optimisationGrid.iterrows():
#        
#        maturity = Utils.monthToYear(row['Maturity'])
#        tenor = Utils.monthToYear(row['Tenor'])
#        
#        # vol is in basis points
#        mktVol = row['Vol']/10000.0
#        
#        # fixed rate is in percentage
#        fixedRate = row['Fss']/100.0
#        print("\n\n")
#        print("Maturity: {}, tenor: {}, vol: {}, fixedRate: {}".format(maturity, tenor, mktVol, fixedRate))
#        
#        normP, normR, G2P, G2R = swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate, mktVol)
#        
#        currErr = (normP + normR - G2P - G2R) ** 2
#        print("Normal Model: Payer: {}, Receiver:{}; G2++ Model: Payer: {}, Receiver: {}".format(normP, normR, G2P, G2R))
#        rmse += currErr
#        print("\n\n")
#    print("Total RMSE: {}".format(rmse))
#    
#g2params = [5.53449534, 3.70068385, 0.001, 0.001, 0.925757890]  
#swaptionModelVsMkt(g2params, optimisationGrid)
    
