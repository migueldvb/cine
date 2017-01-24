# -*- coding: utf-8 -*-

from unittest import TestCase
import cine
import numpy as np

class Test(TestCase):

    def test_is_array(self):
        g = cine.fluorescence.gfactor('H2O')
        self.assertTrue(isinstance(g.gcube, np.ndarray))
