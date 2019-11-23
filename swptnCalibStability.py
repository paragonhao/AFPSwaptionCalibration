#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:50:23 2019

@author: paragonhao
"""


import pandas as pd
from swptnG2PPAF import SWPTNG2PPAF
from utils import Utils
import numpy as np

dates = ['2019-07-05','2019-07-08','2019-07-09','2019-07-10','2019-07-11'
         ,'2019-07-12','2019-07-15','2019-07-16','2019-07-17','2019-07-18','2019-07-19']

optimisationGrid = Utils.readOptimGrid('data/optimisationGrid.csv')

g2params = [2.81905739, 0.00206957183, 0.0479059218, 0.0166433728, 0.998501817]

for date in dates:
    
    marketVolData = pd.read_excel('data/swaption_vols.xlsx',index_col=0,  usecols = "A:L", sheetname=date)
    
    # Todo: get forward starting swap rate 
    fssGrid = pd.read_excel('data/swaption_fss.xlsx',index_col=0,  usecols = "A:L", sheetname=date)
    
    currOptim = optimisationGrid
    currOptim['Vol'] = 0.0
    currOptim['Fss'] = 0.0
    rmse = 0 
    
    for idx, row in optimisationGrid.iterrows():
        
        mktVol = marketVolData[row['Tenor']][row['Maturity']]
        
        currOptim['Vol'][idx] = mktVol
        currOptim['Fss'][idx] =  fssGrid[row['Tenor']][row['Maturity']]

        maturity = Utils.monthToYear(row['Maturity'])
        tenor = Utils.monthToYear(row['Tenor'])
    
#    print("Current Date: {}".format(date))
#    print("Current Optim Grid")
#    print(currOptim)
    
    payFreq = 0.25
    termStructure = Utils.getTermStructure(date)
#    print("Current termStructure at 10 yr {}".format(termStructure.discount(10)))
    
    G2PP = SWPTNG2PPAF(termStructure, payFreq)
    G2PP.swaptionG2PPOptim(g2params, currOptim)

