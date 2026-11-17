"""Setup selection is independent of cuts; metadata never rejects a pulse."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from detector import (CutDatabase, PadGeometry, SETUP_MAPS, channel_to_pad,
                      resolve_pad_map, raw_histogram_spec, raw_2d_axes, load_pad_records,
                      pad_number, pad_rect_mm, scan_run_contexts, default_veto_file, AUXILIARY_DIR, module_floods)


class SetupTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def database(self, text):
        path = self.folder / 'cuts.ini'
        path.write_text(text)
        return CutDatabase(path)

    def test_run_metadata_inherits_and_unknown_is_not_zero(self):
        db = self.database('[global]\nsetup = large_medium\ndrift_voltage_v = 8000\ntilt_angle_deg = 3.5\n'
                           '[run 01616]\nsetup = medium_small\ndrift_voltage_v = -7500\ntilt_angle_deg =\n')
        self.assertEqual(db.run_info(1615), dict(run=1615, setup='large_medium', drift_voltage_v=8000, tilt_angle_deg=3.5))
        self.assertEqual(db.run_info(1616), dict(run=1616, setup='medium_small', drift_voltage_v=-7500, tilt_angle_deg=None))
        self.assertEqual(db.run_info(1616, 'large_medium')['setup'], 'large_medium')
        np.testing.assert_array_equal(db.keep({'Time': np.array([1, 999])}, 1616), [True, True])

    def test_metadata_errors_include_db_line(self):
        for value in ('setup = wrong', 'drift_voltage_v = nan', 'tilt_angle_deg = inf',
                      'drift_voltage_v = 8_000', 'tilt_angle = 0', 'setup = large_medium\nsetup = medium_small'):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, r'cuts.ini:\d+:'):
                self.database('[global]\n'+value+'\n')

    def test_voltage_scan_uses_db_and_available_files_without_rounding(self):
        for name in ('ntuple1.root','ntuple2.root','ntuple3.root','ntuple4.root','ntuple_noise.root'):
            (self.folder / name).touch()
        db = self.database('[global]\ndrift_voltage_v = -8000.5\n'
                           '[run 2]\ndrift_voltage_v =\n'
                           '[run 3]\ndrift_voltage_v = 0\n'
                           '[run 4]\ndrift_voltage_v = 9000.25\n'
                           '[run 99]\ndrift_voltage_v = 1234\n')
        runs, metadata, skipped = scan_run_contexts(self.folder, db, setup='medium_small')
        self.assertEqual(set(runs), {1,3,4})  # 99 has no input; 2 explicitly unknown.
        self.assertEqual([runs[r]['info']['drift_voltage_v'] for r in (1,3,4)], [-8000.5,0,9000.25])
        self.assertTrue(all(r['info']['setup']=='medium_small' for r in runs.values()))
        self.assertEqual(metadata['2']['drift_voltage_v'], None)
        self.assertEqual([r['run'] for r in skipped], [2])
        self.assertFalse((self.folder/'datafiles-voltage.txt').exists())

    def test_skipped_voltage_run_keeps_veto_fallback_but_clean_processed_run_does_not(self):
        live = self.folder/'voltage/rejected_events.csv'
        live.parent.mkdir()
        live.write_text('run,event\n')
        live.with_name('run_summary.csv').write_text('run,noise_events\n1616,0\n')
        with patch('detector.results_directory', return_value=self.folder):
            self.assertEqual(default_veto_file(run=1616), live)
            self.assertEqual(default_veto_file(run=1605), AUXILIARY_DIR/'rejected_events.csv')
            alternate = self.folder/'setups/medium_small/voltage/rejected_events.csv'
            alternate.parent.mkdir(parents=True)
            alternate.write_text('run,event\n2000,10\n')
            alternate.with_name('run_summary.csv').write_text('run\n2000\n')
            self.assertEqual(default_veto_file(run=2000, setup='medium_small'), alternate)
            self.assertEqual(default_veto_file(run=1616, setup='medium_small'), live)

    def test_setup_and_explicit_csv_precedence_reach_pad_cuts(self):
        # Synthetic wiring deliberately swaps two channels between setups.
        # No assumptions about the pending small-pad electronics transcription.
        files = []
        for index, channels in enumerate(((504, 152), (152, 504))):
            path = self.folder / f'map{index}.csv'
            with path.open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['Pad','FADC','Channel','PadSize','Row','Col','X_mm','Y_mm','Width_mm','Height_mm'])
                for pad, channel in enumerate(channels):
                    writer.writerow([pad,channel+1,channel,'Medium',0,pad,4*pad,0,4,4])
            files.append(path)
        db = self.database('[global]\nsetup = large_medium\nreject = Pad == 0\n[run 2]\nsetup = medium_small\n')
        with patch.dict(SETUP_MAPS, dict(large_medium=str(files[0]), medium_small=str(files[1]))):
            rows = {'Channel': np.array([504, 152]), 'Pad': np.array([0, 1])}
            np.testing.assert_array_equal(db.keep(rows, 1), [False, True])
            # Stale supplied Pad values must not override the selected mapping.
            np.testing.assert_array_equal(db.keep(rows, 2), [True, False])
            self.assertEqual(resolve_pad_map(db, 2, 'large_medium'), files[0])
            self.assertEqual(resolve_pad_map(db, 2, 'medium_small', files[0]), files[0])
            np.testing.assert_array_equal(db.keep(rows, 2, files[0]), [False, True])

    def test_pad_axes_use_whole_selected_setup_with_explicit_overrides(self):
        edges, _ = raw_histogram_spec(np.array([-1, 12]), 'Pad', pad_count=508)
        self.assertEqual((len(edges)-1, edges[0], edges[-1]), (508, -.5, 507.5))
        x, y = raw_2d_axes(dict(Pad=np.array([0]), Hit=np.array([1])),
                          dict(X='Pad', Y='Hit', XMax='100', XBins='25'), pad_count=508)
        self.assertEqual((len(x)-1, x[0], x[-1]), (25, -.5, 100))
        edges, _ = raw_histogram_spec(np.array([1]), 'Pad', limits=(0, 50), bins=10, pad_count=508)
        self.assertEqual((len(edges)-1, edges[0], edges[-1]), (10, 0, 50))

    def test_small_grid_blanks_do_not_collapse_positions_or_gain_pad_ids(self):
        # Synthetic channels test geometry independently of electronics approval.
        records = [dict(r) for r in load_pad_records() if r['PadSize'] == 'Medium']
        for pad, record in enumerate(records):
            record.update(Pad=pad, Y_mm=record['Y_mm']-40)
        for row in range(10):
            for col in range(12, 38):
                if (row, col) in {(2, 25), (7, 25)}:
                    continue
                channel = 1000 + len(records)
                records.append(dict(Pad=len(records), FADC=channel+1, Channel=channel,
                                    PadSize='Small', Row=row, Col=col,
                                    X_mm=2*col, Y_mm=40+2*row, Width_mm=2, Height_mm=2))
        path = self.folder / 'synthetic_medium_small.csv'
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader(); writer.writerows(records)
        geometry = PadGeometry(path)
        self.assertEqual((len(geometry), len(geometry.blanks), geometry.width, geometry.height), (508, 242, 100, 60))
        self.assertEqual(pad_number('Small', 0, 12, path), 250)
        self.assertEqual(pad_number('Small', 2, 25, path), -1)
        self.assertEqual(pad_rect_mm('Medium', 0, 0, path), (0, 0, 4, 4))
        np.testing.assert_array_equal(channel_to_pad([r['Channel'] for r in records], path), np.arange(508))
        coverage = np.zeros((60, 100), dtype=int)
        bounds = [tuple(r[k] for k in ('X_mm','Y_mm','Width_mm','Height_mm')) for r in geometry.records] + geometry.blanks
        for x,y,w,h in bounds:
            x,y,w,h = map(int, (x,y,w,h))
            coverage[y:y+h, x:x+w] += 1
        np.testing.assert_array_equal(coverage, np.ones((60,100)))

    def test_installed_medium_small_matches_connected_sketch_and_preserves_other_setup(self):
        path = resolve_pad_map(setup='medium_small')
        geometry = PadGeometry(path)
        self.assertEqual((len(geometry), len(geometry.blanks), geometry.width, geometry.height), (506,244,100,60))
        sections = [r['PadSize'] for r in geometry.records]
        self.assertEqual(sections, ['Medium']*250 + ['Small']*256)
        for fadc, pad in [(153,0),(128,249),(437,250),(290,505)]:
            self.assertEqual(int(channel_to_pad(fadc-1,path)),pad)
        self.assertEqual(int(channel_to_pad(152)),160)  # Existing LM map is unchanged.
        self.assertEqual(int(channel_to_pad(504)),0)
        self.assertEqual(pad_number('Small',2,25,path),-1)
        self.assertEqual(pad_number('Small',3,25,path),-1)
        np.testing.assert_array_equal(channel_to_pad([576,578],path),[-1,-1])
        medium={r['FADC'] for r in geometry.records if r['PadSize']=='Medium'}
        small={r['FADC'] for r in geometry.records if r['PadSize']=='Small'}
        self.assertFalse(medium & small)
        self.assertGreaterEqual(min(small),288)
        self.assertEqual(max(small),575)
        expected=sorted((r['Row'],r['Col']) for r in geometry.records[250:])
        self.assertEqual([(r['Row'],r['Col']) for r in geometry.records[250:]],expected)

    def test_replacement_fadc_picture_audit_matches_every_connected_position(self):
        folder = resolve_pad_map(setup='medium_small').parent
        audit_path = folder/'small_image_audit.csv'
        with audit_path.open(newline='') as stream:
            audit = list(csv.DictReader(stream))
        expected = {(int(r['Row']),int(r['Col'])):int(r['ImageFADC'])
                    for r in audit if r['Connected']=='1'}
        actual = {(r['Row'],r['Col']):r['FADC']
                  for r in load_pad_records(folder/'pad_map_medium_small.csv') if r['PadSize']=='Small'}
        self.assertEqual(len(expected),256)
        self.assertEqual(actual,expected)
        for r in audit:
            if r['Connected']=='0':
                self.assertEqual(r['Pad'],'-1')
                self.assertNotIn((int(r['Row']),int(r['Col'])),actual)
        provenance = json.loads((folder/'pad_map_medium_small_provenance.json').read_text())
        self.assertEqual(provenance['expert_confirmation']['disconnected_fadc'],[577,579])
        exceptions=[r for r in audit if r['ImageFADC'] in ('577','579')]
        self.assertEqual(len(exceptions),2)
        self.assertTrue(all(r['Connected']=='0' and r['Pad']=='-1' for r in exceptions))
        check = provenance['replacement_image_check']
        self.assertEqual(check['connected_positions_checked'],len(expected))
        self.assertEqual(check['audit_sha256'],hashlib.sha256(audit_path.read_bytes()).hexdigest())
        # Source pictures are intentionally not required in portable bundles.
        source = folder.parents[1]/check['source']
        if source.exists():
            self.assertEqual(check['source_sha256'],hashlib.sha256(source.read_bytes()).hexdigest())

    def test_high_channel_module_cannot_merge_with_next_event(self):
        event = np.array([10,10,11,11])
        channel = np.array([576,578,0,1])
        pairs = np.unique(np.column_stack((event,channel)),axis=0)
        self.assertEqual(module_floods(event,channel,np.full(4,50),pairs,3,3),{})
        event = np.array([10,10,10,11,11])
        channel = np.array([576,577,578,0,1])
        pairs = np.unique(np.column_stack((event,channel)),axis=0)
        bad = module_floods(event,channel,np.full(5,50),pairs,3,3)
        self.assertEqual(set(bad),{10})
        self.assertEqual((bad[10][0]['module'],bad[10][0]['channels']),(8,3))

    def test_channel_axes_preserve_requested_spans_and_user_limits(self):
        edges,_ = raw_histogram_spec(np.array([0,152]),'Channel',channel_count=579)
        self.assertEqual((len(edges)-1,edges[0],edges[-1]),(579,-.5,578.5))
        x,y = raw_2d_axes(dict(Channel=np.array([0]),Time=np.array([10])),
                          dict(X='Channel',Y='Time'),channel_count=579)
        self.assertEqual((len(x)-1,x[-1]),(579,578.5))
        edges,_ = raw_histogram_spec(np.array([576,578]),'Channel',channel_count=579,
                                      limits=(-.5,575.5),bins=576)
        self.assertEqual((len(edges)-1,edges[-1]),(576,575.5))

    def test_expert_gap_has_complete_old_new_pad_correspondence(self):
        folder = resolve_pad_map(setup='medium_small').parent
        with (folder/'pad_renumbering_medium_small.csv').open(newline='') as stream:
            migration = list(csv.DictReader(stream))
        self.assertEqual(sorted(int(r['OldPad']) for r in migration),list(range(508)))
        self.assertEqual(sorted(int(r['Pad']) for r in migration if r['Pad']!=''),list(range(506)))
        self.assertEqual({int(r['FADC']) for r in migration if r['Pad']==''},{577,579})
        for r in migration:
            self.assertEqual(int(channel_to_pad(int(r['Channel']),folder/'pad_map_medium_small.csv')),
                             int(r['Pad']) if r['Pad']!='' else -1)


if __name__ == '__main__':
    unittest.main()
