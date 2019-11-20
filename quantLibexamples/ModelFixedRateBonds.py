#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 23:15:19 2019

@author: paragonhao
"""
import QuantLib as ql


todaysDate = ql.Date(5, 7, 2019)

ql.Settings.instance().evaluationDate = todaysDate
spotDates = [ql.Date(15, 1, 2015), ql.Date(15, 7, 2015), ql.Date(15, 1, 2016)]
spotRates = [0.0, 0.005, 0.007]

dayCount = ql.Thirty360()
calendar = ql.UnitedStates()
interpolation = ql.Linear()
compounding = ql.Compounded
compoundingFrequency = ql.Annual
spotCurve = ql.ZeroCurve(spotDates, spotRates, dayCount, calendar, interpolation,
                             compounding, compoundingFrequency)

# construct the term structure curve
spotCurveHandle = ql.YieldTermStructureHandle(spotCurve)

today = spotCurve.referenceDate()
end = today + ql.Period(1, ql.Years)

dates = [ ql.Date(serial) for serial in range(today.serialNumber(), end.serialNumber()+1) ]
spotCurve.enableExtrapolation()

rates_c = [spotCurve.forwardRate(d, ql.TARGET().advance(d,1,ql.Days), ql.Actual360(), ql.Simple).rate() for d in dates ]



issueDate = ql.Date(15, 1, 2015)
maturityDate = ql.Date(15, 1, 2016)
tenor = ql.Period(ql.Semiannual)
calendar = ql.UnitedStates()
bussinessConvention = ql.Unadjusted
dateGeneration = ql.DateGeneration.Backward
monthEnd = False
schedule = ql.Schedule (issueDate, maturityDate, tenor, calendar, bussinessConvention,
                            bussinessConvention , dateGeneration, monthEnd)

list(schedule)

dayCount = ql.Thirty360()
couponRate = .06
coupons = [couponRate]
settlementDays = 0
faceValue = 100
fixedRateBond = ql.FixedRateBond(settlementDays, faceValue, schedule, coupons, dayCount)


# create a bond engine with the term structure as input;
# set the bond to use this bond engine
bondEngine = ql.DiscountingBondEngine(spotCurveHandle)
fixedRateBond.setPricingEngine(bondEngine)

fixedRateBond.NPV()
