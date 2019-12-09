#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 14:18:52 2019

@author: paragonhao
"""


import pandas as pd
from swptnG2PPAF import SWPTNG2PPAF
from utils import Utils
import numpy as np
import matplotlib.pyplot as plt

dates = ['2019-07-05']

file_Path_swaption_vol = 'data/swaption_vols.xlsx'
file_Path_swaption_fss = 'data/swaption_fss.xlsx'
file_Path_swaption_ts = 'data/swaption_termStructure.xlsx'

optimisationGrid = Utils.readOptimGrid(file_Path_optim_grid)
Maturity= ['1M','3M','6M','9M','12M','24M','36M','60M','84M','120M','180M','240M','360M']
Tenor= ['12M','24M','36M','60M','84M','120M','180M','240M','300M','360M']

optimisationGrid = [(x,y) for x in Maturity for y in Tenor]

optimisationGrid = pd.DataFrame(optimisationGrid, columns = ['Maturity', 'Tenor']) 

def reportRMSEOverTime(g2params, dates, reportType, percentage=False):
    rmse_time_series = []

    for date in dates:
        
        marketVolData = pd.read_excel(file_Path_swaption_vol,index_col=0,  usecols = "A:L", sheetname=date)
        
        # Todo: get forward starting swap rate 
        fssGrid = pd.read_excel(file_Path_swaption_fss,index_col=0,  usecols = "A:L", sheetname=date)
        
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
        termStructure = Utils.getTermStructure(date, file_Path_swaption_ts)
    #    print("Current termStructure at 10 yr {}".format(termStructure.discount(10)))
        
        G2PP = SWPTNG2PPAF(termStructure, payFreq)
        
        # True is report type if price
        if reportType:
            if percentage:
                print("RMSE PRICE PERCENTAGE")
                rmse_time_series.append(G2PP.swaptionG2PPOptimPrice(g2params, currOptim, percentage = True))
            else:
                print("RMSE PRICE")
                rmse_time_series.append(G2PP.swaptionG2PPOptimPrice(g2params, currOptim))
        else:
            if percentage:
                print("RMSE VOL PERCENTAGE")
                rmse_time_series.append(G2PP.swaptionG2PPOptimVol(g2params, currOptim, percentage = True))
            else:
                print("RMSE VOL")
                rmse_time_series.append(G2PP.swaptionG2PPOptimVol(g2params, currOptim)) 
    return rmse_time_series


########################################### execution Price RMSE ########################################################
g2params = [2.81905739, 0.00206957183, 0.0479059218, 0.0166433728, 0.998501817]

rmse_over_time = reportRMSEOverTime(g2params, dates, reportType=True, percentage = False)

pd.DataFrame(rmse_over_time).to_csv("rmse_over_time.csv")
rmse_price = pd.read_csv("rmse_over_time.csv")
rmse_price = rmse_price[['0']]


plt.plot(dates, rmse_price)
plt.xticks(rotation=90)
plt.title("Market Price to Model Price RMSE Over 10 days")
plt.ylabel("RMSE")
plt.show()
#############################################################################################################

########################################### execution Volatility ########################################################
g2params = [2.81869666,0.03444719,0.05798171,0.01852287,0.99899906]

rmse_over_time_vol = reportRMSEOverTime(g2params, dates, reportType=False, percentage = False)

pd.DataFrame(rmse_over_time_vol).to_csv("rmse_over_time_vol.csv")
rmse_over_time_vol = pd.read_csv("rmse_over_time_vol.csv")

dates = ['2019-07-05','2019-07-08','2019-07-09','2019-07-10','2019-07-11'
         ,'2019-07-12','2019-07-15','2019-07-16','2019-07-17','2019-07-18','2019-07-19']
rmse_over_time_vol = rmse_over_time_vol[['0']]
plt.plot(dates, rmse_over_time_vol)
plt.xticks(rotation=90)
plt.title("Market Vol to Model Vol RMSE Over 10 days")
plt.ylabel("RMSE")
plt.show()
#############################################################################################################


########################################### execution Price Percentage RMSE ######################################################
g2params = [2.26839563,0.0130107027,0.00106873214,0.0183266914,0.980094659]

rmse_price_percentage = reportRMSEOverTime(g2params, dates, reportType=True, percentage = True)

pd.DataFrame(rmse_price_percentage).to_csv("rmse_price_percentage.csv")
rmse_price_percentage = pd.read_csv("rmse_price_percentage.csv")
rmse_price_percentage = rmse_price_percentage[['0']]
plt.plot(dates, rmse_price_percentage)
plt.xticks(rotation=90)
plt.title("Market Price to Model Price Percentage RMSE Over 10 days")
plt.ylabel("RMSE")
plt.show()
#############################################################################################################

########################################### execution Volatility Percentage ########################################################
g2params = [2.00974568,0.668163150,0.00111277528,0.001,0.988178844]
dates = ['2019-07-05']

rmse_vol_percent = reportRMSEOverTime(g2params, dates, reportType=False, percentage = True)

pd.DataFrame(rmse_vol_percent).to_csv("rmse_vol_percent.csv")
rmse_vol_percent = pd.read_csv("rmse_vol_percent.csv")


rmse_vol_percent = rmse_vol_percent[['0']]
plt.plot(dates, rmse_vol_percent)
plt.xticks(rotation=90)
plt.title("Market Vol to Model Vol RMSE Percentage Over 10 days")
plt.ylabel("RMSE")
plt.show()
#############################################################################################################












