#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:36:52 2019

@author: paragonhao
"""
from g2PlusPlusParams import G2PlusPlusParams as G2Parameters
from g2PlusPlusAnalyticalFormula import G2PlusPlusAnalyticalFormula as GA
import unittest
from mock import patch


class TestStringMethods(unittest.TestCase):
    
    
    def test_g2params_creation(self):
        g2params = G2Parameters(1,1.5,2,0.16,0.5)
        self.assertEqual(g2params.alpha, 1)
        self.assertEqual(g2params.beta, 1.5)
        self.assertEqual(g2params.eta, 2)
        self.assertEqual(g2params.sigma, 0.16)
        self.assertEqual(g2params.rho, 0.5)
     
    def test_ZCB_Var(self):  
        g2params = G2Parameters(1,2,3,4,5)
        var = GA.ZCB_Var(g2params, 0, 1)
        self.assertEqual(round(var,5), 18.68325)
        
    def test_ZCB_Func_A(self):
        g2params = G2Parameters(1,2,3,4,5)
        val = GA.ZCB_Func_A(g2params, 0, 1)
        self.assertEqual(val, 1)
        
if __name__ == '__main__':
    unittest.main()

