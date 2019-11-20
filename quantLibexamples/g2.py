#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 23:14:05 2019

@author: paragonhao
"""

import QuantLib as ql
from collections import namedtuple
import math
from pandas import DataFrame


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

spotCurve.allowsExtrapolation
dates, rates = zip(*spotCurve.nodes())
forwarCurve = ql.ForwardCurve(dates, rates, ql.Actual360())

print(forwarCurve.referenceDate(), 'to', forwarCurve.maxDate())
term_structure = ql.YieldTermStructureHandle(forwarCurve)
index = ql.Euribor3M(term_structure)


# tenor,  maturity, volatility

CalibrationData = namedtuple("CalibrationData", 
                             "start, length, volatility")
#data_normal = [CalibrationData(12, 12, 0.00711),
#        CalibrationData(24, 24, 0.00772),
#        CalibrationData(36, 36, 0.00762),
#        CalibrationData(60, 12, 0.00713),
#        CalibrationData(60, 24, 0.00759),
#        CalibrationData(60, 36, 0.00770),
#        CalibrationData(60, 60, 0.00777),
#        CalibrationData(60, 84, 0.00761),
#        CalibrationData(60, 120, 0.00725),
#        CalibrationData(120, 12, 0.00658),
#        CalibrationData(120, 24, 0.00747),
#        CalibrationData(120, 36, 0.00751),
#        CalibrationData(120, 60, 0.00747),
#        CalibrationData(120, 84, 0.00730),
#        CalibrationData(120, 120, 0.00719),]

data_black = [CalibrationData(12, 12, 0.00711),
        CalibrationData(24, 24, 0.0772),
        CalibrationData(36, 36, 0.0762),
        CalibrationData(60, 12, 0.0713),
        CalibrationData(60, 24, 0.0759),
        CalibrationData(60, 36, 0.0770),
        CalibrationData(60, 60, 0.0777),
        CalibrationData(60, 84, 0.0761),
        CalibrationData(60, 120, 0.0725),
        CalibrationData(120, 12, 0.0658),
        CalibrationData(120, 24, 0.0747),
        CalibrationData(120, 36, 0.0751),
        CalibrationData(120, 60, 0.0747),
        CalibrationData(120, 84, 0.0730),
        CalibrationData(120, 120, 0.0719),]


def create_swaption_helpers_black(data, index, term_structure, engine):
    swaptions = []
    fixed_leg_tenor = ql.Period(3, ql.Months)
    fixed_leg_daycounter = ql.Actual360()
    floating_leg_daycounter = ql.Actual360()
    for d in data:
        vol_handle = ql.QuoteHandle(ql.SimpleQuote(d.volatility))
        helper = ql.SwaptionHelper(ql.Period(d.start, ql.Months),
            ql.Period(d.length, ql.Months),
            vol_handle,
            index,
            fixed_leg_tenor,
            fixed_leg_daycounter,
            floating_leg_daycounter,
            term_structure
        )
        helper.setPricingEngine(engine)
        swaptions.append(helper)
    return swaptions



def create_swaption_helpers_normal(data, index, term_structure, engine):
    
    swaptions = []
    
    fixed_leg_tenor = ql.Period(3, ql.Months)
    
    fixed_leg_daycounter = ql.Actual360()
    
    floating_leg_daycounter = ql.Actual360()
    
    for d in data:
        vol_handle = ql.QuoteHandle(ql.SimpleQuote(d.volatility))
        helper= ql.SwaptionHelper(ql.Period(d.start, ql.Months),
        ql.Period(d.length, ql.Months),
        vol_handle,
        index,
        fixed_leg_tenor,
        fixed_leg_daycounter,
        floating_leg_daycounter,
        term_structure,
        ql.BlackCalibrationHelper.RelativePriceError,
        ql.nullDouble(),
        1.0,
        ql.Normal
    )
    helper.setPricingEngine(engine)
    swaptions.append(helper)
    return swaptions   



def calibration_report(swaptions, data):
    columns = ["Model Price", "Market Price", "Implied Vol", "Market Vol", "Rel Error Price", "Rel Error Vols"]
    report_data = []
    cum_err = 0.0
    cum_err2 = 0.0
    for i, s in enumerate(swaptions):
        
        model_price = s.modelValue()
        
        market_vol = data[i].volatility
        
        black_price = s.blackPrice(market_vol)
        
        rel_error = model_price/black_price - 1.0
        
        implied_vol = s.impliedVolatility(model_price, 1e-5, 50, 0.0, 0.50)
        
        rel_error2 = implied_vol/market_vol-1.0
        
        cum_err += rel_error*rel_error
        
        cum_err2 += rel_error2*rel_error2
        
        report_data.append((model_price, black_price, implied_vol,
                            
        market_vol, rel_error, rel_error2))
        
        print("Cumulative Error Price: %7.5f" % math.sqrt(cum_err))
        print("Cumulative Error Vols : %7.5f" % math.sqrt(cum_err2))
        
    return DataFrame(report_data,columns= columns, index=['']*len(report_data))


model = ql.HullWhite(term_structure);

engine = ql.JamshidianSwaptionEngine(model)

swaptions = create_swaption_helpers_black(data_black, index, term_structure, engine)

optimization_method = ql.LevenbergMarquardt(1.0e-8,1.0e-8,1.0e-8)

end_criteria = ql.EndCriteria(10000, 100, 1e-6, 1e-8, 1e-8)

model.calibrate(swaptions, optimization_method, end_criteria)

a, sigma = model.params()
print("a = %6.5f, sigma = %6.5f" % (a, sigma))
calibration_report(swaptions, data_black)

#
#model = ql.G2(term_structure)
#engine = ql.TreeSwaptionEngine(model, 25)
## engine = ql.G2SwaptionEngine(model, 10, 400)
## engine = ql.FdG2SwaptionEngine(model)
#swaptions = create_swaption_helpers_black(data, index, term_structure, engine)
#
#optimization_method = ql.LevenbergMarquardt(1.0e-8,1.0e-8,1.0e-8)
#end_criteria = ql.EndCriteria(1000, 100, 1e-6, 1e-8, 1e-8)
#model.calibrate(swaptions, optimization_method, end_criteria)
#calibration_report(swaptions, data)
#a, sigma, b, eta, rho = model.params()
#print("a = %6.5f, sigma = %6.5f, b = %6.5f, eta = %6.5f, rho = %6.5f " % (a, sigma, b, eta, rho))