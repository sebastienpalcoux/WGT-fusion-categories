"""Accounting tests for the rank-parameterized driver; no census is run."""
import copy
import sys
import tempfile
from unittest.mock import patch
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run as driver


class DriverAccountingTests(unittest.TestCase):
    def setUp(self):
        self.previous_rank = driver.RANK
        driver.RANK = 44
        self.spec = {'rank': 43, 'sha256': 'a'*64,
                     'tasks': [{'q': 3, 'B': 2, 'm': 1}]}
        self.receipt = {
            'index': 0, 'manifest_sha256': 'a'*64,
            'nodes': 17, 'native': 0, 'handoffs': 0, 'exact_rows': 0,
            'profiles': 0, 'seconds': 0.01, 'native_exit': 0,
            'exact_continuations_exhausted': True,
            'survivors': [], 'excluded': {},
            'raw_ledger': 'rows/000000.jsonl.gz',
            'raw_records_sha256': 'b'*64,
            'native_completion': 'rank 43, q 3, B 2, multiplicity 1, nodes 17, candidates 0, exact 0',
        }

    def tearDown(self):
        driver.RANK = self.previous_rank

    def test_receipt_uses_its_own_rank(self):
        self.assertEqual(driver.check_receipt(self.receipt, self.spec), 0)

    def test_wrong_rank_in_native_trailer_rejected(self):
        rec = copy.deepcopy(self.receipt)
        rec['native_completion'] = rec['native_completion'].replace('rank 43', 'rank 44')
        with self.assertRaises(ValueError):
            driver.check_receipt(rec, self.spec)

    def test_unexhausted_handoff_is_never_complete(self):
        rec = copy.deepcopy(self.receipt)
        rec['exact_continuations_exhausted'] = False
        with self.assertRaises(ValueError):
            driver.check_receipt(rec, self.spec)

    def test_unaccounted_exclusion_rejected(self):
        for reason in ('prime_square_trace', 'prime_cut'):
            rec = copy.deepcopy(self.receipt)
            rec['excluded'] = {reason: 1}
            with self.assertRaises(ValueError):
                driver.check_receipt(rec, self.spec)

    def test_failed_native_exit_is_never_complete(self):
        rec = copy.deepcopy(self.receipt)
        rec['native_exit'] = 2
        with self.assertRaises(ValueError):
            driver.check_receipt(rec, self.spec)

    def test_exact_rows_without_handoffs_rejected(self):
        rec = copy.deepcopy(self.receipt)
        rec.update(exact_rows=1, profiles=1, survivors=[[1,6]])
        with self.assertRaises(ValueError):
            driver.check_receipt(rec, self.spec)

    def test_native_record_count_mismatch_rejected(self):
        rec = copy.deepcopy(self.receipt)
        rec['native'] = 1
        with self.assertRaises(ValueError):
            driver.check_receipt(rec, self.spec)

    def run_incomplete_native(self, program, continuation=None):
        # The stub is only a protocol fixture. No mathematical row is asserted.
        with tempfile.TemporaryDirectory(prefix='rank-accounting-') as temporary:
            work = Path(temporary)
            (work/'rows').mkdir()
            task = {'q':1, 'B':4, 'prefix':[4]}
            with patch.object(driver, 'command', return_value=[sys.executable,'-c',program]):
                with patch.object(driver, 'reason_for', return_value=None):
                    if continuation is None:
                        with self.assertRaises((ValueError, RuntimeError)):
                            driver.run_task(work, work, 0, task)
                    else:
                        with patch.object(driver.square, 'exact_rows', side_effect=continuation):
                            with self.assertRaises((ValueError, RuntimeError)):
                                driver.run_task(work, work, 0, task)
            self.assertFalse((work/'rows/000000.jsonl.gz').exists())
            self.assertFalse(list((work/'rows').glob('*.tmp')))

    def test_truncated_native_record_cannot_publish_a_ledger(self):
        self.run_incomplete_native("import sys; sys.stdout.write('[]')")

    def test_native_failure_after_matching_trailer_cannot_publish(self):
        self.run_incomplete_native(
            "import sys; sys.stderr.write('DONE nodes=1 profiles=0 exact=0'); sys.exit(2)")

    def test_interrupted_exact_generator_cannot_publish_a_ledger(self):
        def incomplete(*args):
            yield [1,6]
            raise RuntimeError('interrupted arbitrary-precision continuation')
        self.run_incomplete_native(
            "import sys; print('{\"exact\":[4]}'); sys.stderr.write('DONE nodes=1 profiles=0 exact=1')",
            incomplete)


    def test_prime_cut_is_declared_by_current_terminal(self):
        import terminal
        self.assertIn('prime_cut', terminal.OBSTRUCTIONS)


if __name__ == '__main__':
    unittest.main()
