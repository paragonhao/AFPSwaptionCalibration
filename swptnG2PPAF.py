#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:51:31 2019

@author: paragonhao
"""
import math 
from utils import Utils
from scipy.stats import norm
import numpy as np
import scipy.optimize as optimize
from normalModel import NormalModel

# G2++ zero coupon bond analytical formula
class SWPTNG2PPAF:
    
    termStructure = None
    notional = 1
    payFreq = None
    
    
    def __init__(self, termStructure, payFreq):
        self.termStructure = termStructure
        self.payFreq = 0.25
        
    
    def getTermStructure(self, T):
        return self.termStructure.discount(T)
    
    
    # V(t,T) in the book page 145
    # g2params: the five parameters 
    # t: starting time of the ZCB 
    # T: ending time of the ZCB
    def V(self, g2params, t, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        Tminust = T - t
        
        expat = math.exp(-alpha * Tminust)
        expbt = math.exp(-beta * Tminust)
        cx = sigma/alpha
        cy = eta/beta
        
        
        
        temp1 = cx * cx * (Tminust + (2.0 * expat - 0.5 * expat * expat - 1.5)/alpha)
        temp2 = cy * cy * (Tminust + (2.0 * expbt - 0.5 * expbt * expbt - 1.5)/beta)
        temp3 = 2.0 * rho * cx * cy * (Tminust + (expat -1.0)/alpha + (expbt - 1.0)/beta - 
                                       (expat * expbt - 1.0)/(alpha + beta))
        
        return temp1 + temp2 + temp3
        
    
    
    # A(t, T) function in the book page 148 
    # g2params: the five parameters 
    # Note: t_i > T
    def A(self, g2params, t, T):
        return (self.getTermStructure(T)/self.getTermStructure(t)) * math.exp(0.5 * (self.V(g2params, t, T) - 
                                      self.V(g2params, 0, T) + self.V(g2params, 0, t)))


    
    def B(self, z, t, T):
        return (1.0 - math.exp(-z*(T - t)))/z
        
    
     # validate from line 69 to 135
    # mu_x function in book page 154
    # mu 1 
    def Mu_x(self, g2params, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        temp = sigma * sigma/(alpha * alpha)
        
        mux = -(temp + rho * sigma * eta/(alpha * beta)) * (1.0 - math.exp(-alpha * T)) \
        + 0.5 * temp * (1.0 - math.exp(-2.0 * alpha * T )) + (rho * sigma * eta/(beta * (alpha + beta))) \
        * (1 - math.exp(-(beta + alpha) * T))
        
        return mux 



    # mu_y function in book page 154
    # mu 2
    def Mu_y(self, g2params, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        temp = eta * eta /(beta * beta)
        
        muy = -(temp + rho * sigma * eta/(alpha * beta)) * (1.0 - math.exp(-beta * T)) \
        + 0.5 * temp * (1.0 - math.exp(-2.0 * beta * T)) \
        + (rho * sigma * eta /(alpha * (alpha + beta))) * (1.0 - math.exp(-(beta + alpha) * T))
        
        return muy
    



    def sigma_x(self, sigma, alpha, T):
        sigma_x = sigma * math.sqrt(0.5 * (1.0 - math.exp(-2.0 * alpha * T))/alpha)
        return sigma_x
    
    
    
    def sigma_y(self, eta, beta, T):
        sigma_y = eta * math.sqrt(0.5 * (1.0 - math.exp(-2.0 * beta * T))/beta)
        return sigma_y
    
    

    def rho_xy(self, g2params, sigma_x, sigma_y, T):
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        return (rho * eta * sigma * (1.0 - math.exp(-(alpha + beta) * T))) /((alpha + beta) * sigma_x * sigma_y)
    
     
################################## Optimisation Function  ############################################
    def swaptionG2PPOptim(self, g2params, optimisationGrid):
        
        rmse = 0.0
             
        print("parameters +++++++++++++++++++++++++++")
        print(g2params)
        print("+++++++++++++++++++++")
        
        for idx, row in optimisationGrid.iterrows():
            
            maturity = Utils.monthToYear(row['Maturity'])
            tenor = Utils.monthToYear(row['Tenor'])
            
            # vol is in basis points
            mktVol = row['Vol']/10000.0

            
            # fixed rate is in percentage
            fixedRate = row['Fss']/100.0
            
            # find out the cash flow and the payment time in terms of years
            # tau: payment Time in terms of year
            # c_i: cash flow at each time of the IRS
            # t_i: time at which the payments are made
            paymentTimes = tenor/self.payFreq
            
            tau = np.repeat(self.payFreq, paymentTimes, axis=0)
            
            c_i = tau * fixedRate
            
            c_i[-1] = c_i[-1] + 1
            t_i = np.arange(maturity+self.payFreq, tenor+maturity + self.payFreq, self.payFreq)
            print("Maturity: {}, tenor: {}, vol: {}, fixedRate: {}".format(maturity, tenor, mktVol, fixedRate))
             
            normP, normR, G2P, G2R = self.swaptionPricingG2PP(g2params, tenor, maturity, self.notional, fixedRate, mktVol, c_i, t_i)
            
            print("normP: {}".format(normP))
            print("normR: {}".format(normR))
            print("G2P: {}".format(G2P))
            print("G2R: {}".format(G2R))
            
            currErr = (normP + normR - G2P - G2R) ** 2
            rmse += currErr
            
        print("Total RMSE: {}".format(rmse))
        return rmse/len(optimisationGrid)
######################################################################################################


######################################################################################################        
    def swaptionPricingG2PP(self, g2params, tenor, maturity, notional, fixedRate, mktVol, c_i, t_i):
    #    print("running the swaptionPricingFunction")
        isPayer = 1
        isReceiver = -1 
        
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        P_0_T = self.getTermStructure(maturity)
        
        sigma_x = self.sigma_x(sigma, alpha, maturity)
        
        sigma_y = self.sigma_y(eta, beta, maturity)
        
        mu_x = self.Mu_x(g2params, maturity)
        
        mu_y = self.Mu_y(g2params, maturity)
        
        rho_xy = self.rho_xy(g2params, sigma_x, sigma_y, maturity)
    
        txy = math.sqrt(1 - rho_xy * rho_xy)
           
        A_array = []
        Ba_array = []
        Bb_array = []
        lambda_array = np.repeat(0, len(t_i), axis=0)
        # define the boundary 
        interval = 0.001
        lower = mu_x - 10 * sigma_x
        upper = mu_x + 10 * sigma_x
        
        xList = np.arange(lower, upper + interval, interval)
        
        for i in range(len(t_i)): 
            
            A_array.append(self.A(g2params, maturity, t_i[i]))
    
            Ba_array.append(self.B(alpha, maturity, t_i[i]))
    
            Bb_array.append(self.B(beta, maturity, t_i[i]))
        
        
        integral_result_Payer = 0 
        integral_result_Receiver = 0 
        cdf_val_1 = 0               # payer
        cdf_val_2 = 0               # receive
        
        # numerically solves the integral
        for x in xList:
            lambda_array = []
            
            for i in range(len(t_i)):
                try:
                    lambda_array.append(c_i[i] * A_array[i] * math.exp(-Ba_array[i] * x))
                except OverflowError:
                    print("exponential value on Ba_array is too big {}")     
            
    
            # confirm if this is correct 
            def y_bar_solver(y):
                sum_all = 0.0
                
                for i in range(len(t_i)):
                    try:
                        sum_all += lambda_array[i] * math.exp(-Bb_array[i] * y)
                    except OverflowError:
                        print("error")
                        return 100
                    
                return sum_all - 1.0
            
            y_bar = optimize.fsolve(y_bar_solver, x0=0, xtol=1e-6)
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
            integral_result_Payer += interval * math.exp(-0.5 * temp_pdf * temp_pdf) * cdf_val_1/(sigma_x * math.sqrt(2.0 * math.pi))
            integral_result_Receiver += interval * math.exp(-0.5 * temp_pdf * temp_pdf) * cdf_val_2/(sigma_x * math.sqrt(2.0 * math.pi))
        
        # swaption price calculated using G2++ analytical solution
        swaptionPricePayer = notional * isPayer * P_0_T * integral_result_Payer
        swaptionPriceReceiver = notional * isReceiver * P_0_T * integral_result_Receiver
        
        # get the discount factor for each payment of the interest rate swap    
        dfs = [self.getTermStructure(t) for t in t_i]
        
        normModel = NormalModel(fRate = fixedRate, strike = fixedRate, impliedVol = mktVol/2.0, maturity = maturity, discountFactors=dfs)
        
        return normModel.getPayerSwaptionPrice(), normModel.getReceiverSwaptionPrice(), swaptionPricePayer, swaptionPriceReceiver