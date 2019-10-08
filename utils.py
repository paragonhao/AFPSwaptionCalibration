#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:06:35 2019

@author: paragonhao
"""
import pandas as pd
import datetime
import math
import numpy as np

class Utils:
    
    @staticmethod
    def readDF(filePath):
        # read in Discount factor
        DFactors = pd.read_csv(filePath)
        DFactors['Date'] = pd.to_datetime(DFactors.Date)
        return DFactors
    
    @staticmethod
    def getZCB_0_T(startDate, maturity, DFactors):
        startDate = datetime.datetime.strptime(startDate, "%m/%d/%y")
        numDays = math.ceil(maturity * 365)
        # find out the time and discount factor 
        P_0_T = DFactors[DFactors['Term[day]'] < numDays].iloc[-1]['DF.mid']
        return P_0_T
    
    @staticmethod
    def getDF(t, DFactors):
        numDays = math.ceil(t * 365)
        return DFactors[DFactors['Term[day]'] < numDays].iloc[-1]['DF.mid']
    
    @staticmethod
    def getDFForT_is(DFactors, tenor, maturity, payFreq):
        t_i = pd.DataFrame(np.arange(maturity+payFreq, tenor+maturity + payFreq, payFreq))
        t_i_df = t_i.apply(lambda x: Utils.getDF(x, DFactors), axis= 1)
        return t_i_df
    
    @staticmethod
    def y_bar_solver(c_i, g2params, t_i_df, P_0_T, t_i, tenor):
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        sum_all = 0 
        for i in range(len(t_i)):
            sum_all += c_i[i] * SWPTNG2PPAF.ZCB_Func_A(g2params, t_i_df.iloc[i], P_0_T, t_i.iloc[i], tenor) * \
            math.exp(-x * SWPTNG2PPAF.B1Formula(alpha, tenor, t_i.iloc[i])) * math.exp( -y_bar * SWPTNG2PPAF.B2Formula(beta, tenor, t_i.iloc[i])) 
        
        return sum_all