#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:51:31 2019

@author: paragonhao
"""
import math 

# G2++ zero coupon bond analytical formula
class SWPTNG2PPAF:
    
    
    # V(t,T) in the book page 145
    # g2params: the five parameters 
    # t: starting time of the ZCB 
    # T: ending time of the ZCB
    @staticmethod
    def ZCB_Variance(g2params, t, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        b1 = SWPTNG2PPAF.B1Formula(alpha, t, T)
        b2 = SWPTNG2PPAF.B2Formula(beta, t, T)
        b12 = SWPTNG2PPAF.B12Formula(alpha, beta, t, T) 
        
        temp1 = ((sigma ** 2)/(alpha ** 2)) * (T - t - b1 - (alpha/2) * (b1 ** 2))
        temp2 = ((eta ** 2)/(beta ** 2)) * (T - t - b2 - (beta/2) * (b2 ** 2))
        temp3 = ((2 * sigma * eta * rho)/(alpha * beta)) * (T - t - b1 - b2 + b12)
        
        return temp1 + temp2 + temp3
        
    
    
    @staticmethod
    def B1Formula(alpha, t, T):
        return (1 - math.exp(-alpha * (T - t)))/alpha
        
        
    
    @staticmethod
    def B2Formula(beta, t, T):
        return (1 - math.exp(-beta * (T - t)))/beta
    
    
    
    @staticmethod
    def B12Formula(alpha, beta, t, T):
        return (1 - math.exp(-(alpha + beta) * (T - t)))/(alpha + beta)
    
    
    
    
    # A(t, T) function in the book page 148 
    # g2params: the five parameters 
    # PM_0_ti: discount factor from 0 to t_i
    # PM_0_T: discount factor from 0 to T
    # Note: t_i > T
    @staticmethod
    def ZCB_Func_A(g2params, PM_0_ti, PM_0_T, t_i, T):
        var_T_ti = SWPTNG2PPAF.ZCB_Variance(g2params, T, t_i)
        var_0_ti = SWPTNG2PPAF.ZCB_Variance(g2params, 0, t_i)
        var_0_T = SWPTNG2PPAF.ZCB_Variance(g2params, 0, T)
        
        return (PM_0_ti/PM_0_T) * math.exp(0.5 * (var_T_ti - var_0_ti + var_0_T))
    
    
    
    # mu_x function in book page 154
    # mu 1 
    @staticmethod
    def Mu_x(g2params, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        b1 = SWPTNG2PPAF.B1Formula(alpha, 0, T)
        b12 = SWPTNG2PPAF.B12Formula(alpha, beta, 0, T)
        
        temp1 = (sigma ** 2)/(2 * (alpha ** 2)) * (1 - math.exp(-2 * alpha * T))
        temp2 = (sigma * eta * rho / beta) * b12
        temp3 = (((sigma ** 2)/alpha) + (sigma * eta * rho/beta)) * b1
        
        return temp1 + temp2 - temp3
    
    
    
    
    # mu_y function in book page 154
    # mu 2
    @staticmethod
    def Mu_y(g2params, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        b2 = SWPTNG2PPAF.B2Formula(beta, 0, T)    
        b12 = SWPTNG2PPAF.B12Formula(alpha, beta, 0, T)
        
        temp1 = (eta ** 2)/(2 * (beta ** 2)) * (1 - math.exp(-2 * beta * T))
        temp2 = (sigma * eta * rho / alpha) * b12 
        temp3 = (((eta ** 2)/beta) + (sigma * eta * rho/alpha)) * b2
        
        return temp1 + temp2 - temp3
    


    @staticmethod 
    def sigma_x(sigma, alpha, T):
        sigma_x = sigma * math.sqrt((1 - math.exp(-2 * alpha * T))/(2 * alpha))
        return sigma_x
    
    
    
    @staticmethod
    def sigma_y(eta, beta, T):
        sigma_y = eta * math.sqrt((1 - math.exp(-2 * beta * T))/(2 * beta))
        return sigma_y
    
    
    
    @staticmethod
    def rho_xy(g2params, sigma_x, sigma_y, T):
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        return ((sigma * eta * rho)/(sigma_x * sigma_y)) * SWPTNG2PPAF.B12Formula(alpha, beta, 0, T)
         
    
    