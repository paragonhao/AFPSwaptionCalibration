#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:11:27 2019

@author: paragonhao
"""


helpers = [ SwapRateHelper(QuoteHandle(SimpleQuote(rate/100.0)),
                           Period(*tenor), TARGET(),
                           Annual, Unadjusted,
                           Thirty360(),
                           Euribor6M())
    for tenor, rate in [((2,Years), 0.201),
                        ((3,Years), 0.258),
                        ((5,Years), 0.464),
                        ((10,Years), 1.151),
                        ((15,Years), 1.588)] ]