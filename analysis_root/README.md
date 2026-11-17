# Native C++ / ROOT analysis

See the [illustrated PDF manual](../docs/TPC_manual.pdf) and editable
[LaTeX source](../docs/TPC_manual.tex). Your `set_color_env()` is called
before drawing in `TPC(...)`; it replaces the standard palette selection.

This is the native C++/ROOT implementation. It does not call
Python, and editing it does not change the working Python analysis.
It requires ROOT with C++17 support; it was tested with ROOT 6.40.04.

## Start here

From the workspace directory containing `analysis_root/` and `data/`:

```bash
root -l 'analysis_root/TPC.C()'
```

Choose from the menu: browse events, save an event, reproduce the examples,
run the voltage analysis, show settings, or check the mapping/cuts.
When the macro finishes you return to the ROOT prompt; `.q` exits ROOT.

You can also call the macro directly:

```bash
# Interactive numbered pads, both medium and large sections.
root -l 'analysis_root/TPC.C("event",1616,5584)'

# Interactive 3D bars instead of the 2D pad plane.
root -l 'analysis_root/TPC.C("event3d",1616,5584)'

# Save an event as 2D/3D PNG/PDF, a hit CSV, and a native ROOT file.
root -l -b -q 'analysis_root/TPC.C("save",1616,5584)'

# Save both established examples and a comparison.
root -l -b -q 'analysis_root/TPC.C("examples")'

# Recompute noise vetoes, select the largest run at each voltage, and plot.
root -l -b -q 'analysis_root/TPC.C("voltage")'

# Fast checks that do not need to read raw data.
root -l -b -q 'analysis_root/TPC.C("check")'
```

During event browsing, **Enter** moves forward, **p** moves back, an event
number jumps, **s** saves the current event, and **q** ends navigation.
Vetoed events and events without accepted mapped hits are skipped.

For a compiled macro, replace `TPC.C(...)` with `TPC.C+(...)`. From an
existing ROOT session:

```cpp
.L analysis_root/TPC.C+
TPC("event", 1616, 5584);
TPC("save", 1616, 5584);
TPC("examples");
TPC("info");
```

## Locations and settings

All paths are resolved from the macro location, not the shell directory.
Default input is `data/root/`. Default output is `results/root/`:

```text
results/root/
  events/run1615/           event plots, hit CSV and event ROOT file
  events/run1616/3d/        3D versions
  comparisons/             the two-run comparison and displayed-hit table
  voltage/                 scan plots, numeric tables and ROOT histograms
```

To use other input/output locations, the full signature is:

```cpp
TPC(command, run, event, inputDirectory, outputDirectory, vetoCSV, cutsDB, studyCutsDB, setup);
```

An empty string keeps the default. For example:

```bash
root -l -b -q '/path/to/analysis_root/TPC.C("save",1616,5584,"/path/to/data","/path/to/output")'
```

Edit **`config/run_cuts.ini`** for global/run-specific manual cuts and
**this folder's `settings.csv`** for automatic-noise thresholds and histogram
limits. The manual DB is shared with Python; the noise settings remain
independent. See [the shared guide](../config/README.md).
The independent voltage comparison additionally reads
`config/voltage_study_cuts.ini`. This file retains its original `Channel < 285`
restriction; event displays and raw/cleaned inspection never read it.
The optional `studyCutsDB` argument overrides that study file only.
It is reread on every `TPC(...)` call. Edit `drift_voltage_v` in each run
block of `config/run_cuts.ini` to set run voltages. This DB is the sole
source; missing/blank voltages are skipped. Run1605 is unknown because
its former table listed both 8250 and 8500 V; resolve it in the DB. Repeated
voltages otherwise select the file with the most distinct events before
cuts; run1603 remains excluded by the settings.

Only two source files are needed:

- `TPC.C`: menu, event displays, exports, voltage scan and audit plots.
- `TPCDetector.h`: raw input, mapping, hit selection, and event-noise logic.

The comments inside both files explain their important operations.
`TPC.C` starts with a reading guide and has ten numbered sections. Start
with Section 10 for command dispatch, Section 2 for drawing styles, and
Sections 6–8 for voltage calculations and output. Function comments explain
inputs, units, ROOT ownership, and the distinction between display hits
and counted raw pulses. Edit manual cuts in `../config/run_cuts.ini`.
`auxiliary/pad_map.csv` is the established 410-pad mapping, including the
corrected large-pad labels. Its physical pad bounds drive the ROOT views.

## Analysis conventions

Raw inspection shows Event, Pad, Hit, Time, Charge, and Pedestal. Pad always
uses the selected setup mapping (large/medium 0–409; medium/small 0–505), independently of cuts. Unconnected
channels are Pad=-1 in underflow, explicitly counted; they remain in the
other raw distributions. `pad_map_used.csv` records the wiring used.
An additional standalone Channel canvas shows all electronics channels
0–575 for large/medium or 0–575 for medium/small, including unconnected ones; printed fADC equals Channel+1. The
six-panel summary retains mapped Pad. DB Channel rules still use electronics
values. Blockhead is omitted. The 1D group opens one 3-by-2 summary,
six individual analysis canvases, and the separate Channel canvas. Its eight
views are saved as PNG/PDF and as editable TCanvas objects in
`raw_histograms.root`, alongside seven TH1 histograms.
`TPC("raw:Channel",1616)` opens only the electronics diagnostic, saving in
`branches_Channel/` without replacing the full run's outputs. `raw:Pad`
selects the mapped-pad diagnostic. Use `root -l results/root/raw/run1616/raw_histograms.root`
and `new TBrowser;` to browse the saved canvases. Python retains equivalent
plots and native ROOT objects for future plotting work.

Event displays use `Time < 228`, `Channel != 271`, `Channel != 149`, and
the run-specific noise veto. They include **both** pad sections. Voltage
spectra read their separate `voltage_study_cuts.ini` for the additional
`Channel < 285` restriction. Channel is zero-based
ROOT electronics numbering; the printed fADC label is `Channel + 1`.
Large pads 0–159 are below medium pads 160–409; pad0 is fADC505. Rows run left to right, then upward. `auxiliary/pad_renumbering.csv` records the old/new correspondence; electronics wiring is unchanged.

Noise is evaluated before all hit cuts. A detector-wide event has at
least 288 distinct channels. Module flooding requires at least 58 distinct
channels in a 72-channel module and either at least 50 channels at one
Time or at least 3 hits per distinct channel. These values are editable.

Displays use the first passing hit per electronics channel, matching the
working Python examples. Voltage spectra count every accepted pulse.
An unmapped channel contributes to time-only spectra, but has no pad plot.
Empty pads draw as zero; negative Charge values are clipped by the display
color scale, not removed from the raw arrays.

The event `.root` file contains a `display_hits` TTree and the native
canvases. It records channel, pad, Hit index, Time, Charge, and transverse
pad centers/widths. It contains the displayed first hits, not every pulse.
The scan `.root` file contains selected-run TH1D/TH2D histograms and plot
canvases; CSVs retain all run summaries and rejected-event evidence.

The 2D pad plane uses TH2Poly. The 3D view uses a 1-mm display grid to draw
ROOT LEGO bars; repeated grid cells are not additional measured hits, so
do not integrate the LEGO histogram for physics. Time heights remain in
native units and are not a reconstructed drift coordinate.

For future reconstruction, start with `ReadRun()` and `PassHit()` in the
header to preserve every raw Event/Channel/Hit/Time/Charge row. The mapping
provides pad bounds and centers. Drift calibration, t0, track fitting and
track output are separate future stages, not implemented here.

## Veto files and sharing

After a native voltage analysis, the other commands prefer its
`results/root/voltage/rejected_events.csv`. Otherwise they use the included
`auxiliary/rejected_events.csv` snapshot. A veto always belongs to a
specific `(run,event)` pair. For new data, regenerate it; for an explicit
ledger, supply the final `vetoCSV` argument. `veto_provenance.json` records
the origin of the bundled snapshot.

Share this folder and its sibling `config/` folder, or
the separate `bundles/analysis_root_bundle.zip`. Raw ROOT data and generated plots
are not part of the source bundle. Compiled `.so`, `.d` and `.pcm` files
can be regenerated locally and are not needed for sharing.

## Raw inspection and manual DB

Normal inspection shows both groups together. Select a group with a flag
inside the command string:

```bash
root -l 'analysis_root/TPC.C("raw",1616)'
root -l 'analysis_root/TPC.C("raw --1d",1616)'
root -l 'analysis_root/TPC.C("raw --2d",1616)'
root -l 'analysis_root/TPC.C("raw-clean",1616)'
root -l 'analysis_root/TPC.C("raw2d:Charge_vs_Pad",1616)'
```

`raw` / `raw-clean` read the run once and show both groups (13 windows).
`--both` is an explicit default. The previous `raw2d` / `raw2d-clean`
commands remain 2D-only shortcuts, and `raw:Branch` stays a focused 1D view.
The 2D group provides a 2-by-2 summary plus four individual canvases: Event vs Pad,
Hit vs Pad, Charge vs Pad, Charge vs Time. The first three have Pad on X;
the current last pair has Time on X. Pair axes/names are editable in CSV.
Colours count Raw Hits. Outputs live in the run's
`2d/` or `cleaned/2d/` subfolder; `raw_2d_histograms.root` contains all four
TH2D objects and five canvases. Edit `config/raw_2d_plots.csv` for pairs and
axis overrides. Current settings use 64 event IDs and inherit one-unit bins
for the other axes. See the
[2D guide](../config/README.md#two-dimensional-inspection) for all options.

Run the synthetic boundary/mapping/flow check with
`root -l -b -q 'analysis_root/tests/Raw2DCheck.C()'`.

```bash
root -l 'analysis_root/TPC.C("raw",1616)'
root -l 'analysis_root/TPC.C("raw:Time",1616)'
root -l 'analysis_root/TPC.C("cuts",1616)'
root -l -b -q 'analysis_root/TPC.C("raw-clean",1616)'
```

Raw means ALL hit rows: no global cut, manual run rule, automatic event
veto, or mapping exclusion. Change per-branch axes in
`config/raw_binning.csv`. Raw outputs are in `results/root/raw/runNNNN/`;
cleaned comparisons are in `cleaned/`. Both keep the full raw axis extent.
The final optional TPC argument selects a different complete cut DB.

Run the detailed rule fixtures with:
`root -l -b -q 'analysis_root/tests/RunCutsCheck.C()'`.

## Run setup metadata

See [setup and run information](../config/README.md#detector-setup-and-run-information).
`setup`, `drift_voltage_v` and `tilt_angle_deg` belong in the detector DB,
with global defaults and per-run overrides. Drift voltage supplies voltage-study assignments; tilt is informational.
Uncut plots still select the run mapping. The additional `auxiliary/pad_map_medium_small.csv` has 506 connected pads:
Medium 0–249 below Small 250–505, retaining 244 unconnected small positions.
Existing runs retain their configured large/medium setup until their DB block is changed.
The replacement `small.jpg` supplies fADC labels, checked at all 256 connected positions;
`small_and_medium.jpg` supplies full-board placement. See `auxiliary/small_image_audit.csv`.

Expert confirmation: 577/579 are disconnected and have no Pad IDs. The central
gap is retained. The 506-pad map supersedes the 508-pad version; use
`auxiliary/pad_renumbering_medium_small.csv` to translate prior Pad rules
and regenerate earlier medium/small exports.

The fADC correspondence diagnostic is `auxiliary/small_vs_large_fadc.pdf`;
its CSV lists all 256 connected Small pads and a Large Pad ID where available. Orange
marks 97 labels with no connected Large-pad counterpart. Disconnected 577/579 are shown as a central gap.
The largest Large label is 576 (Large Pad145). This is a wiring diagnostic.
