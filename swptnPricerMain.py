#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:59:00 2019

@author: paragonhao
"""
import scipy.optimize as optimize
from swptnG2PPAF import SWPTNG2PPAF
from utils import Utils
import pandas as pd



################################## initialize Parameters ############################################
file_Path_swaption_vol = 'data/swaption_vols.xlsx'
file_Path_optim_grid = 'data/optimisationGrid.csv'
file_Path_swaption_fss = 'data/swaption_fss.xlsx'
file_Path_swaption_ts = 'data/swaption_termStructure.xlsx'

g2params = [2.8187, 0.035, 0.0579, 0.0091, 0.999]
calibrationDate = '2019-07-05'
payFreq = 0.25
bnds = ((0.001, 5),(0.001, 5),(0.001, 5),(0.001, 5),(-0.999, 0.999))
######################################################################################################



################################## data initilaisation ############################################
mktVol = pd.read_excel(file_Path_swaption_vol,index_col=0,  usecols = "A:L", sheetname=calibrationDate)
optimisationGrid = Utils.readOptimGrid(file_Path_optim_grid)
fssGrid = pd.read_excel(file_Path_swaption_fss,index_col=0,  usecols = "A:L", sheetname=calibrationDate)

optimisationGrid['Vol'] = 0.0

for idx, row in optimisationGrid.iterrows():
    curVol = mktVol[row['Tenor']][row['Maturity']]
    optimisationGrid['Vol'][idx] = curVol

optimisationGrid['Fss'] = 0.0

for idx, row in optimisationGrid.iterrows():
    optimisationGrid['Fss'][idx] =  fssGrid[row['Tenor']][row['Maturity']]
######################################################################################################



########################################### execution ########################################################
# Contruct the term structure and run the calibration
termStructure = Utils.getTermStructure(calibrationDate, file_Path_swaption_ts)
G2 = SWPTNG2PPAF(termStructure, payFreq)
G2.swaptionG2PPOptim(g2params, optimisationGrid)
result = optimize.minimize(G2.swaptionG2PPOptim, g2params, args=(optimisationGrid), bounds=bnds, method='SLSQP')


if result.success:
    fitted_params = result.x
    print("current fitted value")
    print(fitted_params)
else:
    raise ValueError(result.message)
#############################################################################################################