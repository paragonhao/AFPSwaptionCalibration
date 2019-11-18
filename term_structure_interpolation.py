#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 16:35:54 2019

@author: paragonhao
"""
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from datetime import date

termStructure = pd.read_csv('data/term_structure_interpolation.csv')

x = termStructure['Term[day]']

y = termStructure['Discount']


f = interp1d(x,y)

xnew = np.arange(0, 18270, 1)

import matplotlib.pyplot as plt

plt.plot(x, y, 'o', xnew, f(xnew), '-')
plt.legend(['data', 'linear'], loc='best')
plt.show()

# base date as of 5th July 2019
base = datetime.date(2019,7,5)
numdays = 18270
date_list = [base + datetime.timedelta(days=x) for x in range(numdays)]


term_structure = pd.DataFrame(data={
    'term_days' : xnew,
    'days': date_list,
    'discounts' : f(xnew)
})

term_structure.to_csv('term_structure_interploated.csv', index=False)

