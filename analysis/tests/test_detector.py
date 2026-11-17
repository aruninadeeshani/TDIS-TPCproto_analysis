"""Regression checks: wiring, geometry, selection, and run-specific vetoes.

Run from the data folder: python analysis/tpc.py check
These anchors document the current detector/cut choices; update them
deliberately when changing settings.csv or wiring geometry.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest
import tempfile
import csv
import numpy as np
from detector import load_pad_records, channel_to_pad, pad_rect_mm, pad_number
from detector import passes_hit_cuts, select_hit_mask, rejected_events_for_run
from detector import CutDatabase


class PadMappingTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / 'reference.ini'
        path.write_text('[global]\nreject = Time >= 228\nreject = Channel in 149,271\n')
        self.reference_db = CutDatabase(path)

    def test_wiring_anchors_and_section_transition(self):
        anchors = {505: 0, 487: 1, 431: 159, 153: 160, 151: 161,
                   149: 162, 147: 163, 16: 184, 163: 185, 128: 409, 576: 145, 568: 149}
        for fadc, pad in anchors.items():
            self.assertEqual(int(channel_to_pad(fadc - 1)), pad)

    def test_continuity_round_trip_and_disconnected_channels(self):
        records = load_pad_records()
        self.assertEqual([r['Pad'] for r in records], list(range(410)))
        self.assertEqual(len({r['FADC'] for r in records}), 410)
        np.testing.assert_array_equal(channel_to_pad([r['Channel'] for r in records]), np.arange(410))
        # Printed fADC 1,2,3,129,144 are not present in either wiring diagram.
        np.testing.assert_array_equal(channel_to_pad([-1, 0, 1, 2, 128, 143, 999]), -np.ones(7))

    def test_migration_table_and_row_order_match_current_wiring(self):
        path = Path(__file__).resolve().parents[1] / 'auxiliary/pad_renumbering.csv'
        with path.open(newline='') as stream:
            migration = list(csv.DictReader(stream))
        records = load_pad_records()
        self.assertEqual(sorted(int(r['OldPad']) for r in migration), list(range(410)))
        for r in migration:
            old, new = int(r['OldPad']), int(r['Pad'])
            self.assertEqual(new, old + 160 if old < 250 else old - 250)
            current = records[new]
            self.assertEqual((current['FADC'], current['Channel'], current['PadSize']),
                             (int(r['FADC']), int(r['Channel']), r['PadSize']))
        # All rows retain their left-to-right wiring order, without serpentine reversal.
        for size, count, columns, offset in [('Large', 160, 20, 0), ('Medium', 250, 25, 160)]:
            for index in range(count):
                r = records[offset + index]
                self.assertEqual((r['PadSize'], r['Row'], r['Col']), (size, index // columns, index % columns))
        self.assertEqual((records[0]['X_mm'], records[0]['Y_mm']), (0, 0))
        self.assertEqual((records[160]['X_mm'], records[160]['Y_mm']), (0, 40))

    def test_geometry_covers_board_once_with_large_below_medium(self):
        coverage = np.zeros((80, 100), dtype=int)
        for r in load_pad_records():
            x, y, w, h = pad_rect_mm(r['PadSize'], r['Row'], r['Col'])
            self.assertEqual((x, y, w, h), (r['X_mm'], r['Y_mm'], r['Width_mm'], r['Height_mm']))
            self.assertEqual(pad_number(r['PadSize'], r['Row'], r['Col']), r['Pad'])
            if r['PadSize'] == 'Large':
                self.assertLess(r['Pad'], 160)
                self.assertLess(y, 40)
            else:
                self.assertGreaterEqual(r['Pad'], 160)
                self.assertGreaterEqual(y, 40)
            x, y, w, h = map(int, (x, y, w, h))
            coverage[y:y+h, x:x+w] += 1
        np.testing.assert_array_equal(coverage, np.ones((80, 100)))

    def test_cuts_keep_electronics_meaning(self):
        # A fixed reference selection must not constrain the user's live DB.
        passes = lambda channel,time: bool(self.reference_db.keep({'Channel':[channel], 'Time':[time]})[0])
        self.assertTrue(passes(148, 227))  # fADC149 / pad162
        self.assertFalse(passes(149, 227))  # reference excludes Channel149, pad169
        self.assertFalse(passes(148, 228))
        self.assertTrue(passes(285, 100))
        self.assertTrue(passes(335, 100))
        self.assertTrue(passes(575, 227))
        self.assertFalse(passes(271, 100))
        self.assertEqual(int(channel_to_pad(149)), 169)

    def test_raw_hit_mask_keeps_multiplicity_and_distinct_analysis_profiles(self):
        # Repeated pulses on Channel148 are distinct raw rows. The display
        # may reduce them later, but a future fitter must be able to see both.
        events = np.array([1, 1, 1, 2, 1, 1])
        channels = np.array([148, 148, 335, 335, 149, 335])
        times = np.array([30, 45, 50, 60, 70, 228])
        np.testing.assert_array_equal(select_hit_mask(events, channels, times, [2], cuts_db=self.reference_db),
                                      [True, True, True, False, False, False])
        study = self.reference_db.path.parent / 'study.ini'
        study.write_text('[global]\nreject = Channel >= 285\n')
        np.testing.assert_array_equal(select_hit_mask(events, channels, times, [2], voltage=True, cuts_db=self.reference_db, study_cuts_db=study),
                                      [True, True, False, False, False, False])

    def test_veto_event_numbers_are_local_to_each_run(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'veto.csv'
            path.write_text('run,event\n1615,10\n1616,20\n')
            self.assertEqual(rejected_events_for_run(1615, path), {10})
            self.assertEqual(rejected_events_for_run(1616, path), {20})


if __name__ == '__main__':
    unittest.main()
