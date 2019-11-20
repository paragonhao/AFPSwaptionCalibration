#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:51:08 2019

@author: paragonhao
"""

import math
import utils

import QuantLib as ql

today = ql.Date(11, ql.December, 2012)
ql.Settings.instance().evaluationDate = today


helpers = [
        ql.DepositRateHelper(ql.QuoteHandle(ql.SimpleQuote(rate/100)),
                             ql.Period(1,ql.Days), fixingDays,
                             ql.TARGET(), ql.Following,
                             False, ql.Actual360())
        for rate, fixingDays in [(0.04, 0), (0.04, 1), (0.04, 2)]
]

eonia = ql.Eonia()

helpers += [
        ql.OISRateHelper(2, ql.Period(*tenor),
                         ql.QuoteHandle(ql.SimpleQuote(rate/100)), eonia)
        for rate, tenor in [(0.070, (1,ql.Weeks)), (0.069, (2,ql.Weeks)),
                            (0.078, (3,ql.Weeks)), (0.074, (1,ql.Months))]
]
        
helpers += [
        ql.DatedOISRateHelper(start_date, end_date,
                              ql.QuoteHandle(ql.SimpleQuote(rate/100)), eonia)
        for rate, start_date, end_date in [
                ( 0.046, ql.Date(16,ql.January,2013), ql.Date(13,ql.February,2013)),
                ( 0.016, ql.Date(13,ql.February,2013), ql.Date(13,ql.March,2013)),
                (-0.007, ql.Date(13,ql.March,2013), ql.Date(10,ql.April,2013)),
                (-0.013, ql.Date(10,ql.April,2013), ql.Date(8,ql.May,2013)),
                (-0.014, ql.Date(8,ql.May,2013), ql.Date(12,ql.June,2013))]
]
        
        
helpers += [
        ql.OISRateHelper(2, ql.Period(*tenor),
                         ql.QuoteHandle(ql.SimpleQuote(rate/100)), eonia)
        for rate, tenor in [(0.002, (15,ql.Months)), (0.008, (18,ql.Months)),
                            (0.021, (21,ql.Months)), (0.036, (2,ql.Years)),
                            (0.127, (3,ql.Years)), (0.274, (4,ql.Years)),
                            (0.456, (5,ql.Years)), (0.647, (6,ql.Years)),
                            (0.827, (7,ql.Years)), (0.996, (8,ql.Years)),
                            (1.147, (9,ql.Years)), (1.280, (10,ql.Years)),
                            (1.404, (11,ql.Years)), (1.516, (12,ql.Years)),
                            (1.764, (15,ql.Years)), (1.939, (20,ql.Years)),
                            (2.003, (25,ql.Years)), (2.038, (30,ql.Years))]
]
        
eonia_curve_c = ql.PiecewiseLogCubicDiscount(0, ql.TARGET(), helpers, ql.Actual365Fixed())
eonia_curve_c.enableExtrapolation()
today = eonia_curve_c.referenceDate()
end = today + ql.Period(2,ql.Years)
dates = [ ql.Date(serial) for serial in range(today.serialNumber(), end.serialNumber()+1) ]
rates_c = [ eonia_curve_c.forwardRate(d, ql.TARGET().advance(d,1,ql.Days), ql.Actual360(), ql.Simple).rate() for d in dates ]

 _, ax = utils.plot()

utils.highlight_x_axis(ax)
utils.plot_curve(ax, dates, [(rates_c,'-')], format_rates=True)