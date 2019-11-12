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


class Utils:
    
    
    @staticmethod
    def readDF(filePath):
        # read in Discount factor
        DFactors = pd.read_csv(filePath)
        DFactors['Date'] = pd.to_datetime(DFactors.Date)
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