import os
import sys
sys.path.append(os.path.join(os.getcwd(), "code"))

import unittest

import numpy as np
import pandas as pd

from utils import (
    clean_spike_times,
)


class TestCleanSpikeTimes(unittest.TestCase):

    def test_concatenation(self):
        spike_times = pd.Series([
            [np.array([1,2]), np.array([3])], 
            [np.array([4]), np.array([]), np.array([5, 6])]
        ])

        cleaned = clean_spike_times(spike_times)

        self.assertEqual(len(cleaned), 5)
        np.testing.assert_array_equal(cleaned[0], np.array([1, 2]))
        np.testing.assert_array_equal(cleaned[1], np.array([3]))
        np.testing.assert_array_equal(cleaned[2], np.array([4]))
        np.testing.assert_array_equal(cleaned[3], np.array([]))
        np.testing.assert_array_equal(cleaned[4], np.array([5, 6]))

    def test_empty_arrays(self):
        spike_times = pd.Series([
            [np.array([1,2]), np.array([3])], 
            [],
            [np.array([4, 5]), np.array([])],
        ])

        cleaned = clean_spike_times(spike_times)

        self.assertEqual(len(cleaned), 4)
        np.testing.assert_array_equal(cleaned[0], np.array([1, 2]))
        np.testing.assert_array_equal(cleaned[1], np.array([3]))
        np.testing.assert_array_equal(cleaned[2], np.array([4, 5]))
        np.testing.assert_array_equal(cleaned[3], np.array([]))

    def test_single_row_array(self):
        spike_times = pd.Series([
            [np.array([1, 2]), np.array([3])], 
            [np.array([4, 5])],
        ])

        cleaned = clean_spike_times(spike_times)

        self.assertEqual(len(cleaned), 3)
        np.testing.assert_array_equal(cleaned[0], np.array([1, 2]))
        np.testing.assert_array_equal(cleaned[1], np.array([3]))
        np.testing.assert_array_equal(cleaned[2], np.array([4, 5]))
