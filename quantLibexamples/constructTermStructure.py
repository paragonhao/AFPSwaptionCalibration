#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 16:31:07 2019

@author: paragonhao
"""

import QuantLib as ql
from QuantLib import *
import matplotlib.pyplot as plt 

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

r = SimpleQuote(0.01)

riskFreeCurve = FlatForward(0, TARGET(), QuoteHandle(r), Actual360())