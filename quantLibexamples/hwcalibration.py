#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 23:57:56 2019

@author: paragonhao
"""

import QuantLib as ql
from collections import namedtuple
import math

today = ql.Date(15, ql.February, 2002);
settlement= ql.Date(19, ql.February, 2002);
ql.Settings.instance().evaluationDate = today;

term_structure = ql.YieldTermStructureHandle(
    ql.FlatForward(settlement,0.04875825,ql.Actual365Fixed())
    )

index = ql.Euribor1Y(term_structure)