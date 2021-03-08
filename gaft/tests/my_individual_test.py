#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Test case for Individual.
'''

import unittest

from ..components import MyBinaryIndividual

class MyIndividualTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = True

    def test_my_binary_encoding(self):
        ''' Make sure individual can decode and encode binary gene correctly.
        '''
        indv = MyBinaryIndividual(ranges=[(1, 31)], eps=1)
        indv.init(solution=[3])

        # Test binary chromsome.
        ref_chromsome = [0, 0, 0, 1, 1]
        self.assertListEqual(indv.chromsome, ref_chromsome)

        # Test decode.
        self.assertListEqual(indv.decode(), [3])

        indv = MyBinaryIndividual(ranges=[(1, 31),(1, 31)], eps=1)
        indv.init(solution=[8, 1])

        # Test binary chromsome.
        ref_chromsome = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]
        self.assertListEqual(indv.chromsome, ref_chromsome)

        # Test decode.
        self.assertListEqual(indv.decode(), [8, 1])

    def test_descriptors(self):
        ''' Make sure descriptors can check the parameters correctly.
        '''
        self.assertRaises(TypeError, MyBinaryIndividual, ranges=0.1, eps=0.001)
        self.assertRaises(TypeError, MyBinaryIndividual, ranges=[(0, 1)], eps='asdf')
        self.assertRaises(ValueError, MyBinaryIndividual, ranges=[(0, 1)], eps=10.0)
        self.assertRaises(ValueError, MyBinaryIndividual, ranges=[(0, 1)], eps=[1, 2])

    def test_init(self):
        ''' Make sure the individual can be initialized correctly.
        '''
        indv = MyBinaryIndividual(ranges=[(1, 7)], eps=1)

        # Check chromsome initialization.
        indv.init(chromsome=[0, 1, 1])
        
        self.assertListEqual([0, 1, 1], indv.chromsome)
        self.assertListEqual(indv.solution, [3])

        # Check solution initialization.
        indv.init(solution=[3])
        
        self.assertListEqual(indv.solution, [3])
        self.assertListEqual(indv.chromsome, [0, 1, 1])

    def test_clone(self):
        ''' Make sure individual can be cloned correctly.
        '''
        indv = MyBinaryIndividual(ranges=[(1, 7)],
                                eps=1).init(solution=[4])
        indv_clone = indv.clone()

        self.assertListEqual(indv.chromsome, indv_clone.chromsome)
        self.assertAlmostEqual(indv.solution[0], indv_clone.solution[0], places=2)
        self.assertEqual(indv.ranges, indv_clone.ranges)
        self.assertEqual(indv.eps, indv_clone.eps)

    def test_lower_bound_check(self):
        ''' Make sure we can construct legal individual when individual beyond lower bound.
        '''
        indv = MyBinaryIndividual(ranges=[(1, 7)], eps=[1])
        chromsome = [0, 0, 0]
        indv.init(chromsome)
        indv.bound_check()
        self.assertNotEqual([0, 0, 0], indv.chromsome)
        self.assertIn(indv.solution, [[1], [2], [4]])

if '__main__' == __name__:
    suite = unittest.TestLoader().loadTestsFromTestCase(MyIndividualTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

