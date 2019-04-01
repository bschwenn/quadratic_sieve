from matrix_reducer import *
import unittest
import numpy as np

class MatrixTestCase(unittest.TestCase):

    def test_add_mod_2(self):
        r1 = np.array([1,0,1,0,1])
        r2 = np.array([0,1,0,1,1])
        r3 = add_mod_2(r1,r2)
        self.assertListEqual(r3.tolist(), [1, 1, 1, 1, 0])

    def test_mod_2_rep(self):
        r1 = np.matrix([[3], [7], [4], [1], [2]])
        r1 = mod_2_representation(r1)
        self.assertListEqual(r1.ravel().tolist(), [[1,1,0,1,0]])

    def test_gaussian_elim(self):
        # r1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        # r1 = eliminate(r1)
        # self.assertListEqual(r1.ravel().tolist(), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
        r2 = np.array([[1, 0, 0, 1], [0, 1, 0, 1], [1, 1, 0, 0], [0, 0, 0, 1]])
        r2 = eliminate(r2)
        self.assertListEqual(r2.ravel().tolist(), [1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1])