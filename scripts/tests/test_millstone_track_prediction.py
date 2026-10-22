import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from millstone_track_prediction import fft_grid_search, forecast, predict


class IndependentTrackPrediction(unittest.TestCase):
    def test_sparse_fft_search_matches_direct_complex_sums(self):
        indices=np.array([0,1,3,4,7,8,10,12,13,17,18,20])
        clock=8910;u=(indices-7)*clock/1e6
        rng=np.random.default_rng(25);b=rng.normal(size=len(u))+1j*rng.normal(size=len(u))
        x=np.column_stack([np.ones(len(u)),u]);inverse=np.linalg.inv(x.T@x)
        chirps=np.array([-250.,0.,250.])
        best,nfft=fft_grid_search(u,b,inverse,clock,chirps)
        maximum=-np.inf
        for chirp in chirps:
            f=np.fft.fftfreq(nfft,clock/1e6)
            mod=np.exp(-2j*np.pi*(f[:,None]*u+.5*chirp*u*u))
            rhs=np.column_stack([mod@b,mod@(u*b)])
            explained=np.real(np.sum(np.conj(rhs)*(rhs@inverse.T),axis=1))
            maximum=max(maximum,explained.max())
        self.assertAlmostEqual(best[0],maximum,places=10)
        rhs=x.T@(b*np.exp(-2j*np.pi*(best[1]*u+.5*best[2]*u*u)))
        self.assertAlmostEqual(best[0],np.vdot(rhs,inverse@rhs).real,places=10)

    @staticmethod
    def events():
        keys = 1712587630000000+np.r_[0, np.cumsum(np.tile([8910, 17820], 20))]
        held = int(keys[20])
        t = (keys-held)/1e6
        delay = 3700+3*t+2*t*t
        frequency = -1000-250*t+.5*t*t
        centre = (.1+.02*t)*np.exp(2j*np.pi*(12*t-.5*250*t*t))
        amplitude = centre*np.exp(-2j*np.pi*frequency*(350+delay)/1e6)
        return [dict(sample=int(k), delay_samples=float(d), doppler_hz=float(f), amplitude=a, snr_db=30.)
                for k, d, f, a in zip(keys, delay, frequency, amplitude)], held

    def test_training_recovers_geometry_and_phase_at_an_unseen_pulse(self):
        rows, held = self.events()
        result = forecast(rows, held, 350.)
        self.assertTrue(result['supported'])
        self.assertNotIn(held, result['training_samples'])
        self.assertAlmostEqual(result['delay_samples'], 3700., places=8)
        self.assertAlmostEqual(result['doppler_hz'], -1000., places=8)
        self.assertGreater(result['explained_fraction'], .999)
        expected = .1*np.exp(2j*np.pi*1000*4050/1e6)
        self.assertLess(abs(complex(*result['amplitude'])-expected), .001)

    def test_all_held_measurements_are_excluded_before_association(self):
        rows, held = self.events()
        expected = predict(rows, held, 350.)
        changed = [dict(r) for r in rows]
        changed[20].update(amplitude=1000j, delay_samples=99999., doppler_hz=50000., snr_db=-100.)
        self.assertEqual(predict(changed, held, 350.), expected)

    def test_unseen_clock_phase_and_duplicate_components_are_rejected(self):
        rows, held = self.events()
        sparse = [dict(r, sample=100000+i*20000) for i, r in enumerate(rows)]
        self.assertEqual(forecast(sparse, 110000, 350.)['reason'], 'held_epoch_not_on_training_clock')
        duplicate = [*rows, dict(rows[0])]
        self.assertEqual(forecast(duplicate, held, 350.)['reason'], 'multiple_components_in_training_pulse')


if __name__ == '__main__':
    unittest.main()
