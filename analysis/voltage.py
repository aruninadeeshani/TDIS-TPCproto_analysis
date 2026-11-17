"""Voltage-scan workflow, called only by tpc.py voltage.

General detector cleaning comes from config/run_cuts.ini. This independent
study adds config/voltage_study_cuts.ini; settings.csv controls noise
thresholds. Edit this file for binning, normalization, run selection and
audit plots. The separate study DB retains the original Channel limit.
The detector DB's inherited drift_voltage_v is the sole scan assignment;
available ntuples with unknown voltage are skipped.
Every accepted raw hit contributes; first-hit display reduction is absent.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from detector import (scan_run_contexts, channel_to_pad, load_pad_records, resolve_pad_map, DEFAULT_INPUT_DIR,
    TIME_CUT_MAX, VOLTAGE_CUT, NOISE_CHANNELS, MODULE_CHANNELS,
    MODULE_SYNCHRONOUS, MODULE_HITS_PER_CHANNEL, CHANNELS_PER_MODULE,
    READOUT_MODULES, EXCLUDED_RUNS, select_hit_mask, module_floods, results_directory, DEFAULT_CUTS_DB, DEFAULT_VOLTAGE_STUDY_DB, cut_database, load_hit_arrays)

CUT = VOLTAGE_CUT


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def snapshot_pad_map(path, destination, out):
    """Freeze the effective wiring used by BOTH cut DBs and the pad plots."""
    contents = path.read_bytes()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(contents)
    records = load_pad_records(destination)
    sections = []
    for record in records:
        if not sections or sections[-1]['name'] != record['PadSize']:
            sections.append(dict(name=record['PadSize'], pad_range=[record['Pad'], record['Pad']]))
        else:
            sections[-1]['pad_range'][1] = record['Pad']
    return dict(csv=str(path.resolve()), snapshot=str(destination.relative_to(out)),
                sha256=hashlib.sha256(contents).hexdigest(), pad_count=len(records),
                channel_count=max(CHANNELS_PER_MODULE*READOUT_MODULES, max(r['Channel'] for r in records)+1),
                pad_range=[0, len(records) - 1], sections=sections)


def default_output_directory(source, setups):
    """Keep the established LM destination; isolate a wholly MS scan."""
    base = results_directory(source)
    if set(setups) == {'medium_small'}:
        base = base / 'setups' / 'medium_small'
    return base / 'voltage'


def render_run1610_audit(out, event, time, events, occupancy, base, bad, global_bad, bins, counts):
    # This is the historical run-1610 diagnostic, not a generic track fit.
    # Left: mutually exclusive event classes. Right: detector-only versus
    # full noise vetoes, with identical hit cuts. Event10743 is annotated
    # because its synchronous burst was outside the voltage-scan channels.
    previous_keep = base & ~np.isin(event, global_bad)
    previous_hist, _ = np.histogram(time[previous_keep], bins=bins)
    _, inverse = np.unique(event, return_inverse=True)
    hits_after_cuts = np.bincount(inverse, weights=base, minlength=len(events))
    detector_veto = np.isin(events, global_bad)
    special = (events == 10743) & np.isin(events, bad) & ~detector_veto
    module_veto = np.isin(events, bad) & ~detector_veto & ~special
    good = ~np.isin(events, bad)
    assert np.all(good.astype(int) + detector_veto + module_veto + special == 1)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for mask, color, marker, size, label in [
        (good, "#1559D4", ".", 15, "Kept (passes noise checks)"),
        (detector_veto, "#D00000", "x", 45, "Vetoed: detector-wide flooding"),
        (module_veto, "#E07000", "^", 45, "Vetoed: module flooding"),
        (special, "#9A00B5", "*", 130, "Vetoed: synchronous event 10743"),
    ]:
        axes[0].scatter(occupancy[mask], hits_after_cuts[mask], s=size, color=color,
                        marker=marker, label=f"{label} ({mask.sum():,})", zorder=3)
    if special.any():
        i = np.flatnonzero(special)[0]
        axes[0].annotate("Event 10743: 61 channels in module 5\nfire at exactly the same Time.\n0 hits pass the time/channel cuts.",
                     xy=(occupancy[i], hits_after_cuts[i]), xytext=(145, 205),
                     fontsize=8, color="#6E0080",
                     arrowprops=dict(arrowstyle="->", color="#9A00B5"),
                     bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#9A00B5", alpha=0.95))
    axes[0].set(xlabel="Distinct channels / event (before cuts)",
                ylabel="Hits / event passing time/channel cuts", title="Run 1610: final event classification")
    axes[0].legend(loc="upper left", fontsize=8)
    axes[1].stairs(previous_hist, bins, color="#D00000", label="Detector-wide veto only")
    axes[1].stairs(counts, bins, color="#1559D4", label="All noise vetoes (final)")
    axes[1].set(xlabel="Time (native branch units)", ylabel="Accepted hits / Time count",
                title="Run 1610: effect of module-flooding vetoes")
    axes[1].legend(fontsize=8)
    fig.text(0.5, 0.025, "Time/channel cuts applied to both spectra and to the left-hand hit counts: " + CUT,
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.065, 1, 1))
    fig.savefig(out / "run1610_noise_audit.png", dpi=180)
    fig.savefig(out / "run1610_noise_audit.pdf")
    plt.close(fig)


def add_arguments(parser):
    """Attach voltage-specific options to tpc.py voltage; no separate command."""
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--setup", choices=['large_medium', 'medium_small'],
                        help="Override the detector DB setup for every run in this scan")
    parser.add_argument("--pad-map", type=Path, default=None,
                        help="Explicit wiring CSV override; otherwise resolve the map from each run's setup")
    parser.add_argument("--cuts-db", type=Path, default=DEFAULT_CUTS_DB,
                        help="Detector cuts and run metadata; inherited drift_voltage_v assigns scan voltages")
    parser.add_argument("--study-cuts-db", type=Path, default=DEFAULT_VOLTAGE_STUDY_DB,
                        help="Independent voltage-study restrictions, applied in addition to detector cleaning")
    parser.add_argument("--noise-channels", type=int, default=NOISE_CHANNELS,
                        help="Reject events with at least this many distinct channels before any cuts")
    parser.add_argument("--exclude-run", action="append", type=int, default=EXCLUDED_RUNS.copy(),
                        help="Run excluded from voltage overlays (1603 excluded by default)")
    parser.add_argument("--module-channels", type=int, default=MODULE_CHANNELS,
                        help="Minimum distinct channels in a 72-channel module (about 80%%)")
    parser.add_argument("--module-synchronous", type=int, default=MODULE_SYNCHRONOUS,
                        help="Required synchronous channels, unless module averages >=3 hits/channel")


def run(args):
    """Calculate the existing voltage scan, audits, CSV/NPZ and ROOT histograms."""
    import ROOT  # Load only when the voltage command actually reads ROOT data.
    source = args.input_dir.resolve()
    ROOT.gROOT.SetBatch(True)
    database = cut_database(args.cuts_db)
    study_database = cut_database(args.study_cuts_db)
    global CUT
    CUT = database.describe() + ' && ' + study_database.describe() + "; per-run DB rules"

    # The detector DB is the sole authority for geometry AND scan voltages;
    # the independent study DB contributes only additional hit restrictions.
    setup_override = getattr(args, 'setup', None)
    run_contexts, run_metadata, skipped_runs = scan_run_contexts(
        source, database, setup=setup_override, pad_map=args.pad_map)
    if not run_metadata:
        raise ValueError(f"No ntupleNNNN.root input files found in {source}")
    if not run_contexts:
        raise ValueError(f"No available runs have known drift_voltage_v; set run voltages in {database.path}")
    out = args.output_dir or default_output_directory(
        source, [context['info']['setup'] for context in run_contexts.values()])
    out.mkdir(parents=True, exist_ok=True)
    run_mappings = {}
    results, noise_rows, occupancy_data, pad_data = [], [], {}, {}
    bins = np.arange(-0.5, TIME_CUT_MAX + 0.5, 1.0)  # One native Time count per bin.
    for run, context in run_contexts.items():
        path = context['input_path']
        effective_map = out / 'pad_maps' / f'run{run}_pad_map_used.csv'
        mapping = snapshot_pad_map(context['path'], effective_map, out)
        run_mappings[str(run)] = dict(setup=context['info']['setup'], **mapping)
        pad_count = mapping['pad_count']
        pad_bins = np.arange(-0.5, pad_count + 0.5)
        arrays, input_file = load_hit_arrays(path)
        event, channel, time = (arrays[k] for k in ("Event", "Channel", "Time"))
        events, rows_per_event = np.unique(event, return_counts=True)
        # Count distinct channels, so multiple hits on one channel cannot mimic
        # a detector-wide event. All channels and times participate here.
        pairs = np.unique(np.column_stack((event, channel)), axis=0)
        event_ids, occupancy = np.unique(pairs[:, 0], return_counts=True)
        assert np.array_equal(events, event_ids)
        module_details = module_floods(event, channel, time, pairs,
                                       args.module_channels, args.module_synchronous)
        global_bad = event_ids[occupancy >= args.noise_channels]
        bad = np.union1d(global_bad, list(module_details)).astype(event.dtype)
        is_noise = np.isin(event, bad)
        # Apply hit cuts after identifying whole-event noise. Keep this
        # voltage=True flag here; event displays use voltage=False.
        base = select_hit_mask(event, channel, time, voltage=True, run=run, cuts_db=database,
                               study_cuts_db=study_database, arrays=arrays, pad_map=effective_map)
        keep = base & ~is_noise
        counts, _ = np.histogram(time[keep], bins=bins)
        pads = channel_to_pad(channel, effective_map)
        mapped = keep & (pads >= 0)  # Unconnected electronics have Pad=-1.
        pad_counts = np.bincount(pads[mapped], minlength=pad_count)
        time_pad, _, _ = np.histogram2d(pads[mapped], time[mapped],
                                       bins=(pad_bins, bins))
        assert int(pad_counts.sum()) + int(np.count_nonzero(keep & (pads < 0))) == int(keep.sum())
        pad_data[run] = dict(counts=pad_counts, time_pad=time_pad)
        assert counts.sum() == keep.sum(), "Histogram range omitted accepted times"
        assert not np.intersect1d(event[keep], bad).size
        for ev in bad:
            i = np.searchsorted(events, ev)
            noise_rows.append(dict(run=run, event=int(ev), distinct_channels=int(occupancy[i]),
                                   total_hits=int(rows_per_event[i]),
                                   hits_passing_base_cuts=int(np.count_nonzero(base & (event == ev))),
                                   detector_wide=bool(ev in global_bad),
                                   flooded_modules=json.dumps(module_details.get(int(ev), []))))
        result = dict(**context['info'],
                      voltage=context['info']['drift_voltage_v'],
                      pad_map=mapping['csv'], pad_map_snapshot=mapping['snapshot'],
                      pad_map_sha256=mapping['sha256'], pad_count=pad_count,
                      channel_count=max(mapping['channel_count'], int(np.max(channel, initial=-1))+1),
                      total_events=len(events), clean_events=len(events) - len(bad),
                      accepted_events=len(np.unique(event[keep])), total_hits=len(event),
                      base_cut_hits=int(base.sum()), noise_events=len(bad),
                      additional_module_noise_events=len(np.setdiff1d(bad, global_bad)),
                      noise_hits_after_base_cuts=int(np.count_nonzero(base & is_noise)),
                      accepted_hits=int(keep.sum()), selected=False,
                      mapped_accepted_hits=int(mapped.sum()),
                      unmapped_accepted_hits=int(np.count_nonzero(keep & (pads < 0))),
                      rejected_event_ids=bad.tolist(), histogram=counts.tolist())
        results.append(result)
        occupancy_data[run] = occupancy
        if run == 1610:
            render_run1610_audit(out, event, time, events, occupancy, base, bad,
                                global_bad, bins, counts)
        input_file.Close()
        print(f"Run {run}: {len(events)} events, reject {bad.tolist()}, {keep.sum()} accepted hits", flush=True)

    selected = []
    for voltage in sorted({r['voltage'] for r in results}):
        candidates = [r for r in results if r['voltage'] == voltage and r['run'] not in args.exclude_run]
        if not candidates:
            continue
        # 'Most events' means distinct Event values, never number of hit rows.
        winner = max(candidates, key=lambda r: (r['total_events'], r['clean_events'], -r['run']))
        winner['selected'] = True
        selected.append(winner)

    if not selected:
        raise ValueError("No eligible runs remain; check DB drift_voltage_v values and --exclude-run")

    mixed_setups = len({r['setup'] for r in selected}) > 1
    colors = ["#1559D4", "#D00000", "#008525", "#9A00B5", "#111111", "#E07000", "#009C9C", "#8B4513"]
    for normalized in (False, True):
        fig, ax = plt.subplots(figsize=(13, 7.5))
        for i, result in enumerate(selected):
            values = np.array(result['histogram'], dtype=float)
            if normalized and values.sum():
                values /= values.sum()
            label = f"{result['voltage']:,} V | run {result['run']} | {result['total_events']:,} events"
            if mixed_setups:
                label += f" | {result['setup']}"
            ax.stairs(values, bins, color=colors[i % len(colors)], linestyle="-",
                      linewidth=1.8, label=label)
        # Counts preserves exposure differences; unit-area compares shapes,
        # dividing by accepted HIT count, not by event count or live time.
        ax.set(xlim=(0, TIME_CUT_MAX), xlabel="Time (native branch units)",
               ylabel="Fraction of accepted hits / Time count" if normalized else "Accepted hits / Time count",
               title="TPC time spectra — " + ("unit-area comparison" if normalized else "raw hit counts"))
        ax.grid(alpha=0.2)
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9)
        fig.text(0.08, 0.055, CUT, fontsize=10)
        fig.text(0.08, 0.025, "Noise veto: detector-wide or module flooding, before cuts. One run per voltage; most distinct events.", fontsize=9)
        fig.subplots_adjust(left=0.08, right=0.68, bottom=0.15, top=0.91)
        name = "time_overlay_normalized" if normalized else "time_overlay_counts"
        fig.savefig(out / f"{name}.png", dpi=180)
        fig.savefig(out / f"{name}.pdf")
        plt.close(fig)

    # Pad numbers have meaning only within one setup and wiring snapshot.
    # Keep the time overlay shared, but never overlay different pad maps.
    pad_groups = {}
    for i, result in enumerate(selected):
        pad_groups.setdefault((result['setup'], result['pad_map_sha256']), []).append((i, result))
    pad_plot_groups = []
    for (setup, digest), group in pad_groups.items():
        mapping = run_mappings[str(group[0][1]['run'])]
        pad_count = mapping['pad_count']
        pad_bins = np.arange(-0.5, pad_count + 0.5)
        sections = mapping['sections']
        description = ', '.join(f"{s['name'].lower()} {s['pad_range'][0]}–{s['pad_range'][1]}" for s in sections)
        # Single-map scans retain the historical filenames, including LM.
        suffix = ''
        if len(pad_groups) > 1:
            suffix = f'_{setup}'
            if sum(key[0] == setup for key in pad_groups) > 1:
                suffix += f'_{digest[:12]}'
        pad_plot_groups.append(dict(setup=setup, pad_map_sha256=digest,
                                    runs=[r['run'] for _, r in group],
                                    hits_by_pad=f'hits_by_pad{suffix}', time_vs_pad=f'time_vs_pad{suffix}'))
        fig, ax = plt.subplots(figsize=(13, 6))
        for i, result in group:
            ax.stairs(pad_data[result['run']]['counts'], pad_bins,
                      color=colors[i % len(colors)], linewidth=1.5,
                      label=f"{result['voltage']:,} V | run {result['run']}")
        for section in sections[1:]:
            ax.axvline(section['pad_range'][0] - 0.5, color="#777777", linewidth=0.8)
        ax.set(xlim=(-0.5, pad_count - 0.5), xlabel="Pad number (0-based)", ylabel="Accepted hits / pad",
               title=f"TPC hits in physical pad order ({setup}): {description}")
        ax.legend(loc="upper right")
        ax.grid(alpha=0.15)
        fig.text(0.5, 0.025, "Detector + study cuts and noise vetoes apply. Both cut DBs and pad plots use the same per-run map snapshot.", ha="center", fontsize=9)
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        fig.savefig(out / f"hits_by_pad{suffix}.png", dpi=180)
        fig.savefig(out / f"hits_by_pad{suffix}.pdf")
        plt.close(fig)

        fig, axes = plt.subplots(len(group), 1, figsize=(12, 3 * len(group)), squeeze=False)
        for ax, (_, result) in zip(axes.flat, group):
            mesh = ax.pcolormesh(pad_bins, bins,
                                 pad_data[result['run']]['time_pad'].T, cmap="viridis", shading="flat", rasterized=True)
            for section in sections[1:]:
                ax.axvline(section['pad_range'][0] - 0.5, color="white", linewidth=1)
            ax.set(xlabel="Pad number (0-based)", ylabel="Time (native units)",
                   title=f"{result['voltage']:,} V | run {result['run']}")
            fig.colorbar(mesh, ax=ax, label="Accepted hits")
        fig.suptitle(f"Time versus continuous pad number ({setup}): {description}")
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(out / f"time_vs_pad{suffix}.png", dpi=150)
        fig.savefig(out / f"time_vs_pad{suffix}.pdf")
        plt.close(fig)
    write_csv(out / "pad_histograms.csv",
              [dict(run=r['run'], voltage=r['voltage'], pad=p, count=int(n),
                    setup=r['setup'], pad_map_sha256=r['pad_map_sha256'])
               for r in selected for p, n in enumerate(pad_data[r['run']]['counts'])],
              ['run', 'voltage', 'pad', 'count', 'setup', 'pad_map_sha256'])
    np.savez_compressed(out / "pad_time_histograms.npz",
                        **{f"run{r['run']}": pad_data[r['run']]['time_pad'] for r in selected})

    occupancy_rows = (len(results) + 3) // 4
    fig, axes = plt.subplots(occupancy_rows, 4, figsize=(15, 3.4 * occupancy_rows), squeeze=False)
    for ax, result in zip(axes.flat, results):
        occ = occupancy_data[result['run']]
        channel_count = max(result['channel_count'], int(np.max(occ, initial=0)))
        ax.hist(occ, bins=np.arange(0.5, channel_count+1.5), histtype="step", color="#0072B2")
        ax.axvline(args.noise_channels - 0.5, color="#D55E00", linestyle="--")
        ax.set(yscale="log", xlim=(0, channel_count), title=f"Run {result['run']} | vetoed: {result['noise_events']}",
               xlabel="Distinct channels / event", ylabel="Events")
    for ax in list(axes.flat)[len(results):]:
        ax.set_visible(False)
    fig.suptitle("Event occupancy before cuts — line: detector-wide threshold; veto totals include module flooding")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / "noise_occupancy.png", dpi=150)
    fig.savefig(out / "noise_occupancy.pdf")
    plt.close(fig)

    # Machine-readable outputs retain run/event provenance. Reconstruction
    # should start from raw selected hits plus geometry, not these binned
    # spectra. The ledger is read by both Python and ROOT event displays.
    fields = [key for key in results[0] if key not in ('histogram', 'rejected_event_ids')]
    write_csv(out / "run_summary.csv", [{k: r[k] for k in fields} for r in results], fields)
    write_csv(out / "rejected_events.csv", noise_rows,
              ['run', 'event', 'distinct_channels', 'total_hits', 'hits_passing_base_cuts',
               'detector_wide', 'flooded_modules'])
    hist_rows = [dict(run=r['run'], voltage=r['voltage'], time=int(t), count=int(n))
                 for r in selected for t, n in zip((bins[:-1] + bins[1:]) / 2, r['histogram'])]
    write_csv(out / "histograms.csv", hist_rows, ['run', 'voltage', 'time', 'count'])
    database.snapshot(out / "cuts_used.ini")
    study_database.snapshot(out / "voltage_study_cuts_used.ini")
    metadata = dict(base_cuts=database.describe() + ' && ' + study_database.describe(), manual_cuts_db=str(database.path),
                    voltage_study_cuts_db=str(study_database.path),
                    per_run_rules={str(r["run"]): [rule.expression for rule in database.active(r["run"])] for r in results},
                    study_rules={str(r["run"]): [rule.expression for rule in study_database.active(r["run"])] for r in results}, noise_min_distinct_channels=args.noise_channels,
                    noise_detection="Before all cuts: detector-wide occupancy OR module flooding",
                    module_noise_rule=dict(channels_per_module=CHANNELS_PER_MODULE, min_distinct_channels=args.module_channels,
                                           min_synchronous_channels=args.module_synchronous,
                                           alternative_min_hits_per_distinct_channel=MODULE_HITS_PER_CHANNEL,
                                           logic="module coverage AND (synchrony OR repeated hits)"),
                    event_count_definition="Number of distinct Event values present in the ntuple, before cuts",
                    pad_mapping=dict(per_run=run_mappings,
                                     electronics_conversion="ROOT Channel = printed one-based FADC - 1",
                                     unmapped="Pad=-1; no pad bin; retained in time-only spectra if it passes both cut DBs"),
                    setup_override=setup_override,
                    pad_map_override=str(args.pad_map.resolve()) if args.pad_map is not None else None,
                    pad_plot_groups=pad_plot_groups,
                    voltage_assignment_source=str(database.path),
                    voltage_assignment_snapshot='cuts_used.ini',
                    voltage_assignment_definition="Detector DB drift_voltage_v, inherited global -> run; unknown voltages skipped; setup override changes geometry only",
                    run_metadata=run_metadata, skipped_runs=skipped_runs,
                    excluded_runs=args.exclude_run, runs=results)
    (out / "analysis.json").write_text(json.dumps(metadata, indent=2) + "\n")
    root_out = ROOT.TFile(str(out / "time_histograms.root"), "RECREATE")
    ROOT.TNamed('analysis_metadata', json.dumps(metadata)).Write()
    for r in selected:
        hist = ROOT.TH1D(f"time_run{r['run']}", f"{r['voltage']} V, run {r['run']} | {r['setup']};Time (native branch units);Accepted hits", TIME_CUT_MAX, -0.5, TIME_CUT_MAX - 0.5)
        for i, count in enumerate(r['histogram'], start=1):
            hist.SetBinContent(i, count)
            hist.SetBinError(i, np.sqrt(count))
        hist.SetEntries(r['accepted_hits'])
        hist.Write()
        ph = ROOT.TH1D(f"pad_run{r['run']}", f"{r['voltage']} V, run {r['run']} | {r['setup']};Pad number (0-based);Accepted hits", r['pad_count'], -0.5, r['pad_count'] - 0.5)
        for i, count in enumerate(pad_data[r['run']]['counts'], start=1):
            ph.SetBinContent(i, int(count))
            ph.SetBinError(i, np.sqrt(count))
        ph.SetEntries(r['mapped_accepted_hits'])
        ph.Write()
    root_out.Close()
    print("Selected:", [(r['voltage'], r['run']) for r in selected])
    print("Skipped runs (unknown drift_voltage_v):", [r['run'] for r in skipped_runs])
    print("Output:", out)
