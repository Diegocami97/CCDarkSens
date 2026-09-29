# CCDarkSens

I built CCDarkSens to compute dark matter sensitivity projections and exclusion limits for DAMIC-M — upper limits on the interaction strength between dark matter and ordinary matter, derived from a profile-likelihood analysis over what a DAMIC-M-style silicon CCD would actually observe. It's a C++17 framework, with the physics that feeds it (rate calculations) written in Python.

If you're new here, this file is meant to be a map, not a manual. It tells you what's where and which document to open next — the actual explanations live in `docs/`.

---

## What this framework actually does, in one paragraph

Start with a dark matter model (a mass and a coupling strength) and a signal-rate table computed for it. Fold that rate through a model of the detector — diffusion, readout noise, charge quantization, pixel-pattern classification, or a full per-event cluster fit, depending on which channel you're running. Combine the result with an expected background. Run a profile-likelihood-ratio scan over the coupling at fixed mass to find the 90% CL upper limit. Repeat across a mass grid. That's the whole pipeline — the interesting part is how each detector-response and background piece is actually built, which is what the code in `src/` and `include/` does.

Four physics channels are implemented: DM-electron scattering (QEDark/QCDark2), dark-photon absorption, the Migdal effect, and WIMP-nucleon spin-independent scattering. Three observable spaces are supported: `n_e` bins, SRDM pattern bins, and (for the WIMP-nucleon channel) reconstructed per-event cluster energy.

---

## Where things live

```
apps/       C++ executables — scans, plotting, validation, one-off diagnostics
src/        the framework's implementation (config parsing, detector/experiment
            setup, signal models, detector response, backgrounds, statistics)
include/    headers for everything in src/
python/     the actual rate physics (ccdarkphys/) — this is where dR/dE gets
            computed from first principles, before any C++ code sees it
utils/      standalone scripts that call into python/ to generate the rate
            CSVs the C++ side reads
configs/    JSON configs — configs/examples/ is the curated, documented set;
            everything else in configs/ feeds one of the apps above directly
data/       reference inputs (ionization tables, efficiency tables, digitized
            literature curves) — only the files actually needed are tracked
docs/       everything explained below
```

A run is always: write or copy a JSON config → point it at a rate directory → run one of the apps in `apps/` → read the ROOT file it writes. Nothing here needs anything besides a config file and a build.

---

## Building it

```bash
cmake -B build -S .
cmake --build build -j8
```

You need ROOT (built with Minuit2) and `nlohmann_json` on your system; CMake will tell you clearly if either is missing. Everything links against the same core library, so a full build gives you every app in `apps/` at once.

---

## Where to actually go next

I'd read these in this order:

1. **[`docs/Beginners_Guide.md`](docs/Beginners_Guide.md)** — start here. Concepts, every config field, worked examples, a glossary. If you only read one document, read this one.
2. **[`docs/Framework_Architecture.md`](docs/Framework_Architecture.md)** — once you've run something and want to know *how* the code actually does it: config parsing, the detector-response layer, the background layer, the statistics layer, module by module, with the exact mechanics of how a config turns into an exclusion curve.
3. **[`docs/DM_Signal_Models_Physics_Reference.md`](docs/DM_Signal_Models_Physics_Reference.md)** — the physics itself: the DM-electron, dark-photon, and Migdal rate calculations, charge ionization, and how they're validated.
4. **[`docs/ClusterFitMC_Design.md`](docs/ClusterFitMC_Design.md)** — the WIMP-nucleon channel specifically: the per-event cluster-fit reconstruction, its noise-tail calibration, and an honest record of where that channel's reproduction of a published result still has an open discrepancy.
5. **[`configs/examples/README.md`](configs/examples/README.md)** and **`docs/Student_Examples_*.md`** — runnable, documented examples for every channel, each with a real-data reproduction and a hypothetical-exposure projection side by side. This is the fastest way to see the whole pipeline actually run.
6. **[`docs/GenericScanApp_Design.md`](docs/GenericScanApp_Design.md)** — if you're going to modify `ccdarksens_scan_generic` itself (the app every new example should use), read this first; it's the record of how that app was built and validated against the older, frozen reference apps.

The `*_Reproduction_Guide.md` files in `docs/` (QEDark, QCDark2, Migdal, dark photon, WIMP-nucleon) are older, narrower step-by-step guides for specific reproductions — still accurate, just more specific in scope than the student examples above.

---

## A note on scope

This repo only tracks what's needed to build the framework and reproduce the examples documented above. A lot of exploratory work — investigation notes, abandoned plans, one-off pydme crosschecks, scratch config sweeps — happened along the way and stays on my own disk rather than in git. If a script or path mentioned somewhere seems to point at something that isn't here, that's why; it wasn't meant to be part of what gets shared.
