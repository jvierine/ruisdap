import importlib.util
import os
import pathlib
import unittest

import numpy as np

root = pathlib.Path(__file__).resolve().parents[2]
path = root / 'scripts/check-millstone-residual-controls.py'
spec = importlib.util.spec_from_file_location('residual_controls', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PipeControls(unittest.TestCase):
    def test_native_pipe_preserves_unseeded_data_and_removes_injected_echo(self):
        binary = pathlib.Path(os.environ.get('MILLSTONE_CANCEL_HELPER', str(root/'target/debug/examples/cancel_in_memory')))
        self.assertTrue(binary.is_file(), 'Build cargo example cancel_in_memory before this integration test')
        interpolation = module.load('test_interpolation', root/'scripts/check-millstone-residual-cfar.py')
        rng = np.random.default_rng(19)
        y = (rng.standard_normal(1024)+1j*rng.standard_normal(1024))/np.sqrt(2)
        tx = np.zeros(len(y), complex)
        tx[16:80] = np.tile([1., 1., -1., 1.], 16)
        native = module.NativeCancellation(binary.resolve())
        try:
            clean, mask, qa = native.cancel(y, tx, [], [800, 1000])
            np.testing.assert_array_equal(clean, y)
            self.assertTrue(mask.all())
            self.assertEqual(qa['reason'], 'no_detections')
            target = (10+4j)*interpolation.template(tx, 250.25, 2000., 1e6, 16)
            clean, mask, qa = native.cancel(y+target, tx, [dict(delay_samples=250.25, doppler_hz=2000.)], [800, 1000])
            self.assertTrue(qa['accepted'])
            self.assertTrue(qa['applied'])
            self.assertLess(np.sum(abs(clean-y)**2), .01*np.sum(abs(target)**2))
            self.assertTrue(mask.all())
        finally:
            native.close()

    def test_paired_moments_use_no_wraparound_or_masked_samples(self):
        y = np.array([1., 2., 3., 1000.], complex)
        np.testing.assert_allclose(module.moments(y, np.array([1, 1, 1, 0], bool), [0, 1, 2]), [14/3, 4., 3.])


if __name__ == '__main__':
    unittest.main()
