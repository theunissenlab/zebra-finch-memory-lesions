import os
import sys
sys.path.append(os.path.join(os.getcwd(), "code"))

import unittest
import warnings

import numpy as np

from stats import (
    _odds_ratio,
    false_discovery,
    jackknife
)


class TestStats(unittest.TestCase):

    def test_false_discovery(self):
        """Test false_discovery implements the Handale-Anscombe correction"""
        p = np.array([0.05])
        np.testing.assert_array_equal(
            false_discovery(p),
            np.array([True])
        )

        p = np.array([0.04, 0.05])
        np.testing.assert_array_equal(
            false_discovery(p),
            np.array([True, True])
        )

        p = np.array([0.04, 0.06])
        np.testing.assert_array_equal(
            false_discovery(p),
            np.array([False, False])
        )

        p = np.array([0.02, 0.06])
        np.testing.assert_array_equal(
            false_discovery(p),
            np.array([True, False])
        )

        p = np.array([0.06, 0.02])
        np.testing.assert_array_equal(
            false_discovery(p),
            np.array([False, True])
        )

        p = np.array([0.06, 0.02, 0.06])
        np.testing.assert_array_equal(
            false_discovery(p),
            np.array([False, False, False])
        )

    def test_odds_ratio_zero_correction(self):
        a = np.array([
            [5, 4],
            [7, 2]
        ])
        OR, se = _odds_ratio(a, zero_correction=True)
        np.testing.assert_almost_equal(OR, 0.4074074)
        np.testing.assert_almost_equal(se, 0.9681806)

        a = np.array([
            [0, 4],
            [7, 2]
        ])
        OR, se = _odds_ratio(a, zero_correction=True)
        np.testing.assert_almost_equal(OR, 0.03703704)
        np.testing.assert_almost_equal(se, 1.65998661)

        a = np.array([
            [1, 0],
            [7, 0]
        ])
        OR, se = _odds_ratio(a, zero_correction=True)
        np.testing.assert_almost_equal(OR, 0.2)
        np.testing.assert_almost_equal(se, 2.19089023)

    def test_odds_ratio_no_zero_correction(self):
        """Test that _odds_ratio doesn't break when no-zero correction is requested

        Tests for divide by zero runtime warnings since they are expected.
        """
        a = np.array([
            [5, 4],
            [7, 2]
        ])

        with warnings.catch_warnings(record=True) as warnings_:
            OR, se = _odds_ratio(a, zero_correction=False)
            self.assertEqual(len(warnings_), 0)

        np.testing.assert_almost_equal(OR, 0.3571429)
        np.testing.assert_almost_equal(se, 1.0453981)

        a = np.array([
            [0, 4],
            [7, 2]
        ])

        with self.assertWarns(RuntimeWarning):
            OR, se = _odds_ratio(a, zero_correction=False)

        np.testing.assert_almost_equal(OR, 0.0)
        np.testing.assert_almost_equal(se, np.inf)

        a = np.array([
            [1, 0],
            [7, 0]
        ])

        with self.assertWarns(RuntimeWarning):
            OR, se = _odds_ratio(a, zero_correction=False)

        np.testing.assert_almost_equal(OR, np.nan)
        np.testing.assert_almost_equal(se, np.inf)

    def test_jackknife(self):
        np.random.seed(420)

        a = np.random.random(size=50)
        m, se = jackknife(a, np.mean)

        true_mean = np.mean(a)
        true_se = np.std(a) / np.sqrt(len(a))

        np.testing.assert_almost_equal(m, true_mean)
        np.testing.assert_almost_equal(se, true_se)
