import importlib.util
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
