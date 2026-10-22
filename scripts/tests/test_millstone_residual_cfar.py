import importlib.util
import pathlib
import unittest

import numpy as np

path = pathlib.Path(__file__).resolve().parents[1] / 'check-millstone-residual-cfar.py'
spec = importlib.util.spec_from_file_location('residual_cfar', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ResidualReconstruction(unittest.TestCase):
    def test_delay_is_non_circular_and_doppler_has_physical_sign(self):
        tx = np.zeros(64, complex)
        tx[4:8] = 1
        expected = np.zeros(64, complex)
        expected[24:28] = np.exp(2j * np.pi * 10000 * np.arange(24, 28) / 1e6)
        np.testing.assert_allclose(module.template(tx, 20., 10000., 1e6, 16), expected, atol=1e-15)
        np.testing.assert_array_equal(module.template(tx, 100., 0., 1e6, 16), np.zeros(64))

    def test_recorded_amplitude_reproduces_residual_and_detects_tampering(self):
        config = dict(instrument_metadata={'dc_offset': [.2, -.3]},
                      lpi={'tx_gate': [16, 32], 'sample_rate_hz': 1e6},
                      cancellation={'interpolation_half': 16})
        tx = np.zeros(256, complex)
        tx[16:32] = np.exp(1j * np.arange(16))
        v = module.template(tx, 90.25, -2500., 1e6, 16)
        y = tx + (3+2j) * v
        support = abs(v) > abs(v).max() * 1e-3
        rms = np.sqrt(np.mean(abs(tx[support])**2))
        diagnostic = dict(detections=[{'delay_samples': 90.25, 'doppler_hz': -2500.}],
                          amplitudes=[[3., 2.]], noise_rms=1., residual_rms_over_noise=rms)
        clean, value = module.reconstruct(y + (.2-.3j), config, diagnostic)
        np.testing.assert_allclose(clean, tx, atol=1e-14)
        self.assertAlmostEqual(value, rms, places=12)
        diagnostic['amplitudes'][0][0] += .01
        with self.assertRaises(AssertionError):
            module.reconstruct(y + (.2-.3j), config, diagnostic)


if __name__ == '__main__':
    unittest.main()
