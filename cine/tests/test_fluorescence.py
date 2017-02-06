# -*- coding: utf-8 -*-

import unittest
import cine
import numpy as np

class Test(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        """
        Handle the arguments using TestCase
        """
        super(Test, self).__init__(*args, **kwargs)
        self.g = cine.fluorescence.pumping('HDO', 7)

    def test_is_array(self):
        """
        Demonstrates that the gfactor attribute is an array object
        """
        self.assertTrue(isinstance(self.g.gfactor, np.ndarray))

    def test_nlevels(self):
        """
        Demonstrates the gfactor size and shape
        """
        self.assertTrue(len(self.g.gfactor), 2)
        self.assertTrue(self.g.gfactor.shape, (2, 7))

    def test_float_dtype(self):
        """
        Demonstrates the gfactor dtype
        """
        self.assertEqual(self.g.gfactor.dtype, np.float64)

    def test_values(self):
        """
        Demonstrates the values in the gfactor array
        """
        ghdo = np.zeros((7, 7))
        ghdo[0, 3] = 2.568872e-05
        ghdo[0, 4] = 2.570305e-05
        ghdo[0, 5] = 1.552757e-05
        ghdo[1, 2] = 6.253229e-05
        ghdo[1, 6] = 2.987896e-05
        ghdo[2, 1] = 6.196215e-05
        ghdo[2, 6] = 4.410062e-05
        ghdo[3, 0] = 7.547422e-05
        ghdo[3, 4] = 3.103947e-05
        ghdo[3, 5] = 5.048423e-05
        ghdo[4, 0] = 1.253741e-04
        ghdo[4, 3] = 5.128064e-05
        ghdo[4, 5] = 4.679292e-05
        ghdo[5, 0] = 7.481781e-05
        ghdo[5, 3] = 8.287649e-05
        ghdo[5, 4] = 4.643613e-05
        ghdo[6, 1] = 4.820172e-05
        ghdo[6, 2] = 7.201329e-05
        self.assertIsNone(np.testing.assert_almost_equal(self.g.gfactor, ghdo, decimal=5))

if __name__ == '__main__':
    unittest.main()
