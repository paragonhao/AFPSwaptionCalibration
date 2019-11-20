#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 12:34:02 2019

@author: paragonhao
"""

import pandas as pd
import datetime
import math
import numpy as np
import re
import QuantLib as ql
import matplotlib.pyplot as plt 

class Utils:
    
    
    @staticmethod
    def readDF(filePath):
        # read in Discount factor
        DFactors = pd.read_csv(filePath)
        DFactors['days'] = pd.to_datetime(DFactors.days)
        return DFactors
    
    @staticmethod 
    def readMktVolSurface(filePath):
        mktVol = pd.read_csv(filePath,index_col=0)
        return mktVol
    
    @staticmethod
    def readOptimGrid(filePath):
        optimGrid = pd.read_csv(filePath)
        return optimGrid
    
    @staticmethod
    def readFSS(filePath):
        fssGrid = pd.read_csv(filePath,index_col=0)
        return fssGrid
    
    @staticmethod
    def monthToYear(monthString):
        monthInYr = 12.0
        month = re.search("(\d+)", monthString)
        return float(month[0])/monthInYr


    @staticmethod
    def getTermStructure():
        todaysDate = ql.Date(5, 7, 2019)
        
        ql.Settings.instance().evaluationDate = todaysDate
        
        spotDates = [ql.Date(10, 7, 2019), ql.Date(18, 7, 2019), ql.Date(25, 7, 2019),ql.Date(1, 8, 2019)
        , ql.Date(13, 8, 2019), ql.Date(11, 9, 2019), ql.Date(11, 10, 2019), ql.Date(14, 11, 2019), ql.Date(11, 12, 2019)
        , ql.Date(13, 1, 2020), ql.Date(13, 4, 2020), ql.Date(13, 7, 2020), ql.Date(13, 1, 2021), ql.Date(13, 7, 2021)
        , ql.Date(13, 7, 2022), ql.Date(12, 7, 2023), ql.Date(11, 7, 2024), ql.Date(13, 7, 2026), ql.Date(11, 7, 2029)
        , ql.Date(11, 7, 2031), ql.Date(12, 7, 2034), ql.Date(13, 7, 2039), ql.Date(13, 7, 2044), ql.Date(13, 7, 2049)
        , ql.Date(11, 7, 2059), ql.Date(11, 7, 2069)]
        
        spotRates = [2.42, 2.386449337, 2.387049675 ,2.388999462 ,2.322000027 ,2.228999615 ,2.160995007, 2.097999573, 2.06099987, 2.015614986
        , 1.920724392, 1.849999428, 1.733134747, 1.662094593, 1.593999863, 1.574999809,1.578999996,1.638650937,1.745031829, 1.811209614, 1.880565407
        , 1.945579175, 1.967107654, 1.973911374, 1.956311315, 1.923911183]
        
        spotRates = [x/100 for x in spotRates]
        
        dayCount = ql.Thirty360()
        calendar = ql.UnitedStates()
        interpolation = ql.Linear()
        compounding = ql.Compounded
        compoundingFrequency = ql.Annual
        
        spotCurve = ql.ZeroCurve(spotDates, spotRates, dayCount, calendar, interpolation,
                                     compounding, compoundingFrequency)
        return spotCurve
    
    