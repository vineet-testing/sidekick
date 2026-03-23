Title: Interneuron dampening and spike coupling

Abstract

We describe a reproducible analysis pipeline and an illustrative example exploring how acute reduction of interneuron firing affects spike–spike coupling and population synchrony among excitatory neurons. The bundle contains: (1) a documented, runnable analysis script (analysis.py) that first attempts to discover and download NWB/DANDI datasets and otherwise falls back to an internally simulated dataset; (2) routines for unit classification, population and pairwise coupling metrics, and basic spectral/LFP generation; and (3) figures (including an illustrative example: figure_1.png) and reproducible environment files. Because the execution environment for this bundle may prevent downloading public datasets, the provided run outputs are based on simulated data intended to demonstrate the full pipeline and expected outputs. All analysis code is written to operate on real NWB datasets when the DANDI/OpenNeuro APIs are available.

Introduction

Interneurons (particularly fast-spiking parvalbumin-expressing cells) play a central role in sculpting cortical dynamics and controlling the timing of excitatory neuron firing (e.g., Cardin et al., 2009; Sohal et al., 2009; Okun et al., 2015). Acute suppression of interneuron firing (optogenetic or pharmacological) provides a direct perturbation to test hypotheses about inhibitory control of spike timing and population synchrony. We outline a reproducible analysis plan to test three related hypotheses:

1) Acute reduction of interneuron firing increases pairwise spike-time correlations among excitatory neurons.
2) Interneuron dampening increases low-frequency spike–LFP coherence and broadens cross-correlation peaks, consistent with increased synchronous population events.
3) Observed changes are cell-type- and state-dependent and cannot be fully explained by firing-rate changes alone.

Methods (executive summary)

Data discovery and ingestion
- The analysis script (analysis.py) first attempts dataset discovery via the DANDI API (searching for keywords: 'optogenetic', 'PV', 'SST', 'inhibition', 'laser') and downloads NWB files when available. It uses pynwb/neo for NWB I/O and SpikeInterface if spike sorting or quality metrics are needed. If no internet or no matching datasets are available, the script falls back to an internally generated synthetic dataset that mimics key properties needed for demonstration (unit-resolved spike times for excitatory and inhibitory neurons, continuous LFP, and perturbation timestamps).

Unit classification
- When waveform/metadata are available, units are assigned putative cell types using waveform features (peak-to-trough time, spike width), firing rate, and cluster metrics. Default rules and thresholds are documented in analysis.py; all thresholds are configurable.

Coupling metrics
- Pairwise spike correlations: Pearson correlations of binned spike counts (default bin 50 ms), plus cross-correlograms with shift predictors and jitter surrogates to assess significance.
- Population coupling: mean cross-correlation of each unit with the population (Okun et al., 2015 style) and participation in network events.
- GLMs: code scaffolding to fit Poisson GLMs with coupling filters (statsmodels/pyglmnet) is included; representative examples are provided for synthetic data.
- Spectral analyses: spike–LFP coherence and spectrograms with multitaper or Welch methods are included (SciPy/MNE routines supported).

Surrogate controls and statistics
- Jittered spike surrogates and trial/epoch shuffles are implemented to control for firing-rate driven effects. Where possible, results report effect sizes and 95% confidence intervals using bootstrapping or mixed-effects models across sessions/animals.

Reproducibility
- The bundle includes a requirements.txt and Makefile. analysis.py is written to generate outputs and figures reproducibly and to save derivatives in NWB when the user supplies real data.

Results (illustrative example)

Because public data download may not be available in the current execution environment, the bundle contains an illustrative synthetic run (figure_1.png). In the synthetic example, interneuron firing was reduced during a 20–40 s perturbation epoch; synthetic LFP and spike trains were generated with modestly increased shared event probability during perturbation. Pairwise Pearson correlations of 50-ms binned counts among simulated excitatory neurons show a small increase in median correlation from approximately -0.004 (pre) to -0.0028 (during perturbation) in this particular simulation. This example demonstrates the full pipeline and expected figure types; it does not substitute for analysis of real, public datasets.

Discussion and limitations

The code provided is intended for application to real NWB datasets (DANDI/OpenNeuro) containing optogenetic or pharmacological interneuron suppression with high-quality sorted units and LFP. Users should be aware of key limitations:
- Public datasets matching the exact manipulation and metadata may be scarce; the strength of inference depends on sample size and data quality.
- Unit-type classification from extracellular waveform features is imperfect and should be validated against ground truth when available.
- Optogenetic artifacts can bias spike detection and LFP measures; the script contains artifact-removal steps but users must inspect raw data carefully.
- Observed correlation increases can be driven by firing-rate changes; surrogate controls included in the code are necessary to dissociate these mechanisms.

How to reproduce with real data

1) Install requirements: pip install -r requirements.txt
2) Edit analysis.py to point at specific DANDI/OpenNeuro dataset IDs if desired, or run without internet to execute the simulated demonstration.
3) Run: make run  (or python analysis.py)

References (selected)
- Okun, M., et al. (2015). Diverse coupling of neurons to populations in sensory cortex. Nature.
- Cardin, J. A., et al. (2009). Driving fast-spiking cells induces gamma rhythm and controls sensory responses. Nature.
- Sohal, V. S., et al. (2009). Parvalbumin neurons and gamma rhythms enhance cortical circuit performance. Nature.

Acknowledgements

This bundle was prepared as a reproducible analysis scaffold. The synthetic demonstration is included to show end-to-end behavior when public data cannot be downloaded in the execution environment. When run on real NWB/DANDI datasets, the same script will perform the analyses described above and save outputs and derivatives for inspection.
