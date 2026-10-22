import importlib.util
import copy
import pathlib
import subprocess
import sys
import tempfile
import unittest

import digital_rf
import h5py
import numpy as np

path = pathlib.Path(__file__).resolve().parents[1] / 'audit-millstone-source-coverage.py'
spec = importlib.util.spec_from_file_location('source_coverage', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class IndependentCoverage(unittest.TestCase):
    def test_insufficient_calibration_exclusion_replays_every_pulse(self):
        keys=list(range(18));valid={2,7}
        gap=dict(reason='insufficient_receiver_calibration',mode='alternating-code-16',pulse_samples=keys,
                 science_pulse_samples=keys[:-1],calibration_pulse_samples=sorted(valid),
                 configuration=dict(mode='alternating-code-16',cancel_debris=True,injection_temperature_k=1200,
                     instrument_metadata={'dc_offset':[0.,0.]},cancellation={'interpolation_half':16},
                     lpi=dict(sample_rate_hz=1e6,filter_half=100,passband_hz=18000,
                              science_pulse=[True]*17+[False],noise_gate=[7700,8200],
                              injection_gate=[8400,8850],tx_gate=[76,624])),
                 cancellation_diagnostics=[dict(reason='no_detections' if k in valid else 'invalid_background_noise',
                                                detections=[],amplitudes=[]) for k in keys])
        class OriginalRF:
            def read_vector(self,key,length,channel):
                return np.exp(2j*np.pi*10000*np.arange(length)/1e6) if key in valid else np.zeros(length,complex)
        np.testing.assert_array_equal(module.verify_calibration_exclusion(gap,OriginalRF()),keys[:-1])
        changed=copy.deepcopy(gap);changed['calibration_pulse_samples']=[2]
        with self.assertRaisesRegex(AssertionError,'silently omitted'):
            module.verify_calibration_exclusion(changed,OriginalRF())
        changed=copy.deepcopy(gap);changed['calibration_pulse_samples']=[2,7,8]
        with self.assertRaises(AssertionError):module.verify_calibration_exclusion(changed,OriginalRF())

    def test_bad_gain_exclusion_requires_original_negative_calibration(self):
        n=10000;omega=2*np.pi*10000/1e6;half=100
        j=np.arange(2*half);m=j-half+1e-6
        fir=(.5-.5*np.cos(2*np.pi*j/(2*half-1)))*np.sin(2*np.pi*21600/1e6*m)/(np.pi*m);fir/=fir.sum()
        response=abs(np.sum(fir*np.exp(1j*omega*(half-j))))**2
        keys=list(range(18))
        gap=dict(reason='invalid_receiver_gain',mode='alternating-code-16',pulse_samples=keys,
                 science_pulse_samples=keys,calibration_pulse_samples=keys,noise_power=4*response,
                 injection_power=response,receiver_gain=-3*response/1200,injection_temperature_k=1200,
                 configuration=dict(mode='alternating-code-16',cancel_debris=True,injection_temperature_k=1200,
                     instrument_metadata={'dc_offset':[0.,0.]},cancellation={'interpolation_half':16},
                     lpi=dict(sample_rate_hz=1e6,filter_half=100,passband_hz=18000,science_pulse=[True]*len(keys),
                              noise_gate=[7700,8200],injection_gate=[8400,8850],tx_gate=[76,624])),
                 cancellation_diagnostics=[dict(reason='no_detections',detections=[],amplitudes=[]) for _ in keys])
        class OriginalRF:
            def __init__(self,injection=1.):
                self.raw=2*np.exp(1j*omega*np.arange(n));self.raw[8300:9000]*=injection/2
            def read_vector(self,key,length,channel):
                self.assertions=(length,channel);return self.raw.copy()
        np.testing.assert_array_equal(module.verify_calibration_exclusion(gap,OriginalRF()),keys)
        with self.assertRaises(AssertionError):module.verify_calibration_exclusion(gap,OriginalRF(3.))
        changed=copy.deepcopy(gap);changed['calibration_pulse_samples']=keys[1:]
        with self.assertRaisesRegex(AssertionError,'silently omitted'):module.verify_calibration_exclusion(changed,OriginalRF())

    def test_whole_vectors_and_internal_gaps_are_required(self):
        keys = np.array([100, 110, 120, 130, 140, 150], np.uint64)
        np.testing.assert_array_equal(module.covered(keys, {100: 25, 140: 20}, length=10),
                                      [True, True, False, False, True, True])

    def test_real_partners_can_cross_the_integration_edge(self):
        keys = np.arange(0, 180000, 10000, dtype=np.uint64)
        codes = np.repeat(np.arange(9), 2)
        result = module.paired_science(keys, codes, 10000, 170000)
        np.testing.assert_array_equal(result, keys[1:-1])
        self.assertEqual(len(module.paired_science(keys, codes, 10000, 160000)), 0)

    def test_unsupported_or_distant_partners_are_not_invented(self):
        keys = np.arange(20, dtype=np.uint64) * 40000
        self.assertEqual(len(module.paired_science(keys, np.ones(20), 0, 1000000)), 0)

    def test_upstream_metadata_and_indexed_gap_detect_missing_science(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            source = root / 'original'
            for channel in ['zenith-l', 'tx-h']:
                folder = source / 'rf_data' / channel
                folder.mkdir(parents=True)
                with digital_rf.DigitalRFWriter(str(folder), np.dtype('complex64'), 3600, 1000,
                                                 0, 1000000, 1, compression_level=1,
                                                 is_complex=True, is_continuous=False) as writer:
                    writer.rf_write(np.zeros((220000, 1), np.complex64))
                    writer.rf_write(np.zeros((770000, 1), np.complex64), next_sample=230000)
            keys = np.arange(20, dtype=np.uint64) * 20000 + 100000
            metadata = {
                'id_metadata': (1000000, keys, {'sweepid': np.repeat(np.arange(1, 11), 2)}),
                'antenna_control_metadata': (1000000, [0], {'rx_antenna': 'ZENITH', 'tx_antenna': 'ZENITH'}),
                'powermeter': (1, [0], {'zenith_power': 500000.}),
            }
            for name, (rate, times, data) in metadata.items():
                folder = source / 'metadata' / name
                folder.mkdir(parents=True)
                digital_rf.DigitalMetadataWriter(str(folder), 3600, 1, rate, 1, 'md').write(times, data)
            audits = root / 'campaign/1970-01-01/audit'
            audits.mkdir(parents=True)
            audit = audits / '0.h5'
            # The gap drops sample 220000, and its real partner at 240000
            # cannot be substituted with a neighbouring pulse of another code.
            expected = keys[(keys != 220000) & (keys != 240000)]
            with h5py.File(audit, 'w') as h:
                h['candidate_science_pulse_sample'] = keys
                h['science_pulse_sample'] = expected
            command = [sys.executable, str(path), '--root', str(source), '--campaign', str(root/'campaign'),
                       '--date', '1970-01-01', '--seconds', '10', '--report', str(root/'report.h5')]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with h5py.File(root/'report.h5') as h:
                self.assertTrue(h.attrs['verified'])
                self.assertEqual(h.attrs['science_pulse_count'], 18)
                self.assertFalse(h.attrs['complete_utc_day'])
            resumed = subprocess.run([*command, '--resume'], capture_output=True, text=True)
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            with h5py.File(audit, 'a') as h:
                h['science_pulse_sample'][0] += 1
            changed = subprocess.run([*command, '--resume'], capture_output=True, text=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn('Previously verified audit changed', changed.stderr)


if __name__ == '__main__':
    unittest.main()
