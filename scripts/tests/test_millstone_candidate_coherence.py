import importlib.util
import pathlib
import unittest

import numpy as np

folder = pathlib.Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, folder/filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


phase = load('phase_controls', 'check-millstone-candidate-coherence.py')
summary = load('summary_controls', 'summarize-millstone-residual-controls.py')
independent = load('independent_controls', 'check-millstone-independent-amplitude.py')


class ScientificControlSummaries(unittest.TestCase):
    def test_pulse_local_fit_phase_is_referred_to_the_echo_epoch(self):
        rng = np.random.default_rng(103)
        times = np.linspace(0, .2, 20)+rng.uniform(0, .0005, 20)
        centres = np.linspace(3700., 3710., 20)
        epoch = times+centres/1e6
        frequencies = 1300.+250.*epoch+rng.normal(0, 60., len(times))
        # Per-pulse Doppler errors shift the fitted zero-epoch amplitude phase.
        amplitudes = np.exp(2j*np.pi*(1300.*epoch+.5*250.*epoch**2-frequencies*centres/1e6))
        result, null = phase.phase_search(times, amplitudes, frequencies, centres, 16, rng)
        self.assertGreater(result['normalized_coherence'], .9)
        self.assertLessEqual(result['random_phase_p'], 2/17)
        self.assertTrue(np.isfinite(null).all())

    def test_paired_bootstrap_keeps_known_bias_and_rejects_zero_reference(self):
        ref = np.arange(1., 11.)
        np.testing.assert_allclose(summary.paired_ratio(.1*ref, ref, np.random.default_rng(4)), [.1, .1, .1])
        self.assertTrue(np.isnan(summary.paired_ratio(ref, np.zeros(10), np.random.default_rng(4))).all())

    def test_held_pulse_amplitude_cannot_change_its_own_prediction(self):
        coefficients = np.array([[1+2j, 5+6j], [2+3j, 6+7j], [3+4j, 7+8j]])
        expected = independent.leave_one_out(coefficients)
        changed = coefficients.copy()
        changed[1] += 10-20j
        np.testing.assert_allclose(independent.leave_one_out(changed)[1], expected[1])
        np.testing.assert_allclose(expected[1], (coefficients[0]+coefficients[2])/2)


if __name__ == '__main__':
    unittest.main()
