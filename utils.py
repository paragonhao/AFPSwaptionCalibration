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
    def readOptimGrid(filePath):
        optimGrid = pd.read_csv(filePath)
        return optimGrid
    
    
    
    @staticmethod
    def monthToYear(monthString):
        monthInYr = 12.0
        month = re.search("(\d+)", monthString)
        return float(month[0])/monthInYr

    
    
    @staticmethod
    def getTermStructure(date):
        
        tsRaw = pd.read_excel('data/swaption_termStructure.xlsx',index_col=None,  usecols = "A,C", sheetname=date)
        
        todaysDate = ql.DateParser.parseFormatted(date,'%Y-%m-%d')
        
        ql.Settings.instance().evaluationDate = todaysDate
        
        spotDates =[]
        spotRates = []
        
        for idx, row in tsRaw.iterrows():
            spotDates.append(ql.DateParser.parseFormatted(row['Payment Date'],'%m/%d/%Y'))  
            spotRates.append(row['Market Rate'])
        
        spotRates = [x/100 for x in spotRates]
        
        dayCount = ql.Thirty360()
        calendar = ql.UnitedStates()
        interpolation = ql.Linear()
        compounding = ql.Compounded
        compoundingFrequency = ql.Annual
        
        spotCurve = ql.ZeroCurve(spotDates, spotRates, dayCount, calendar, interpolation,
                                     compounding, compoundingFrequency)
        return spotCurve
    