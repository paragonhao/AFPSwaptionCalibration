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
import QuantLib as ql

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
    def getTermStructure(date, filePath):
        
        tsRaw = pd.read_excel(filePath, index_col=None, usecols = "A,C", sheetname=date)
        
        todaysDate = ql.DateParser.parseFormatted(date, '%Y-%m-%d')
        
        ql.Settings.instance().evaluationDate = todaysDate
        
        spotDates =[]
        spotRates = []
        
        for idx, row in tsRaw.iterrows():
            spotDates.append(ql.DateParser.parseFormatted(row['Payment Date'], '%m/%d/%Y'))  
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
    
    
    
    @staticmethod
    def getForwardCurve(spotCurve):
        today = spotCurve.referenceDate()
        end = today + ql.Period(50,ql.Years)
        
        dates = [ ql.Date(serial) for serial in range(today.serialNumber(), end.serialNumber()+1) ]
        
        rates_c = [spotCurve.forwardRate(d, ql.TARGET().advance(d,1,ql.Days), ql.Actual360(), ql.Simple).rate() for d in dates ]
        return rates_c