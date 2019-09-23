#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:51:31 2019

@author: paragonhao
"""
import math 

class G2PlusPlusAnalyticalFormula:
    
    # waiting for DF data to coming in
    # should be in a separate class
    @staticmethod
    def getMarket_DF(T):
        
        return 1
    
    # V(t,T) in the book page 145
    # g2params: the five parameters 
    # t: starting time of the ZCB 
    # T: ending time of the ZCB
    @staticmethod
    def ZCB_Var(g2params, t, T):
        temp_1 = ((g2params.sigma / g2params.alpha) ** 2) \
                    * (T - t + (2/g2params.alpha) * math.exp(-g2params.alpha * (T - t)) \
                    - (1 / (2 * g2params.alpha)) * math.exp(-2 * g2params.alpha * (T - t)) \
                    - (3/(2 * g2params.alpha)))
        
        temp_2 = ((g2params.eta / g2params.beta) ** 2) * \
                    (T - t + (2/g2params.beta) * math.exp(-g2params.beta * (T - t)) \
                    - (1 / (2*g2params.beta)) * math.exp(-2 * g2params.beta * (T - t)) \
                    - (3/(2 * g2params.beta)))
                    
        temp_3 = 2 * g2params.rho * ((g2params.sigma * g2params.eta)/(g2params.alpha * g2params.beta))\
                * (T - t + (math.exp(-g2params.alpha * (T - t)) - 1)/g2params.alpha \
                   + (math.exp(-g2params.beta * (T - t)) - 1)/g2params.beta \
                   - (math.exp(-(g2params.alpha + g2params.beta)*(T - t)) - 1)/(g2params.alpha + g2params.beta))
                
        return temp_1 + temp_2 + temp_3
    
    
    # A(t, T) function in the book page 148 
    # g2params: the five parameters 
    # t: starting time of the ZCB 
    # T: ending time of the ZCB
    @staticmethod
    def ZCB_Func_A(g2params, t, T):
        P_M_0_T = G2PlusPlusAnalyticalFormula.getMarket_DF(T)
        P_M_0_t = G2PlusPlusAnalyticalFormula.getMarket_DF(t) 
        return (P_M_0_T/P_M_0_t) * math.exp(0.5 * (G2PlusPlusAnalyticalFormula.ZCB_Var(g2params, t, T) \
        -  G2PlusPlusAnalyticalFormula.ZCB_Var(g2params, 0, T) +  G2PlusPlusAnalyticalFormula.ZCB_Var(g2params, 0, t)))
    
    # B(z, t, T) function in the book page 148 
    # g2params: the five parameters 
    # z: variable 
    # t: starting time of the ZCB 
    # T: ending time of the ZCB
    @staticmethod
    def ZCB_Func_B(z, t, T):
        temp1 =  1 - math.exp(-z * (T-t))
        return temp1/z
    
    # M(s, t) function for x in the book page 154
    # s: starting time of the ZCB  
    # t: starting time of the option
    # T: ending time of the process
    @staticmethod
    def Moment_T_x(g2params, s, t, T):
        temp1 = ((g2params.sigma/g2params.alpha)**2) + (g2params.rho * g2params.sigma * g2params.eta/(g2params.alpha * g2params.beta))
        temp2 = (1 - math.exp(-g2params.alpha * (t-s))) * temp1
        temp3 = 0.5 * ((g2params.sigma/g2params.alpha)**2) * (math.exp(-g2params.alpha * (T-t)) - math.exp(-g2params.alpha * (T+t-2*s)))        
        temp4 = ((g2params.rho * g2params.sigma * g2params.eta) / (g2params.beta * (g2params.alpha + g2params.beta))) \
        * (math.exp(-g2params.beta * (T-t)) - math.exp(-g2params.beta * T - g2params.alpha * t + (g2params.alpha + g2params.beta) * s))     
        return temp2 - temp3 - temp4
    
    # M(s, t) function for y 
    @staticmethod
    def Moment_T_y(g2params, s, t, T):
        temp1 = ((g2params.eta/g2params.beta)**2) + (g2params.rho * g2params.sigma * g2params.eta/(g2params.alpha * g2params.beta))
        temp2 = (1 - math.exp(-g2params.beta * (t-s))) * temp1
        temp3 = 0.5 * ((g2params.eta/g2params.beta)**2) * (math.exp(-g2params.beta * (T-t)) - math.exp(-g2params.beta * (T+t-2*s)))        
        temp4 = ((g2params.rho * g2params.sigma * g2params.eta) / (g2params.alpha * (g2params.alpha + g2params.beta))) \
        * (math.exp(-g2params.alpha * (T-t)) - math.exp(-g2params.alpha * T - g2params.beta * t + (g2params.alpha + g2params.beta) * s))     
        return temp2 - temp3 - temp4
    
    @staticmethod 
    def sigma_x(g2params, T):
        sigma_x = g2params.sigma * math.sqrt((1 - math.exp(-2 * g2params.alpha * T))/(2 * g2params.alpha))
        return sigma_x
    
    @staticmethod
    def sigma_y(g2params, T):
        sigma_y = g2params.eta * math.sqrt((1 - math.exp(-2 * g2params.beta * T))/(2 * g2params.beta))
        return sigma_y
    
    @staticmethod
    def rho_xy(g2params, T, sigma_x, sigma_y):
        rho_xy = ((g2params.rho * g2params.sigma * g2params.eta)/((g2params.alpha + g2params.beta) * sigma_x * sigma_y)) \
        * (1 - math.exp(-(g2params.alpha + g2params.beta) * T))
        return rho_xy
    
    