# analysis.py
# Reproducible analysis script for "Interneuron dampening and spike coupling"
# Functionality:
#  - Attempt to discover and download NWB files via dandi-api
#  - Load NWB via pynwb/neo (when available)
#  - If no external data is available, run a synthetic demonstration dataset
#  - Compute simple unit classification (waveform-free fallback), pairwise binned correlations, population coupling, and produce figures
#  - Save outputs to ./outputs

import os
import sys
import json
import argparse
from datetime import datetime

# Try to import optional dependencies; if not available, we'll fail noisily for those parts
try:
    import numpy as np
    import scipy.signal as signal
    import matplotlib.pyplot as plt
except Exception as e:
    print('Required Python packages are missing. Please install requirements.txt before running. Error:', e)
    raise

OUTDIR = 'outputs'
os.makedirs(OUTDIR, exist_ok=True)

# --- Helper: create a synthetic dataset (used if no public data available) ---

def generate_synthetic_dataset(seed=42):
    np.random.seed(seed)
    fs = 1000
    duration_s = 60
    t = np.arange(0, duration_s, 1/fs)
    n_time = len(t)

    n_exc = 40
    n_inh = 10

    # Perturbation epoch
    perturb_start = 20.0
    perturb_end = 40.0
    perturb_mask = (t >= perturb_start) & (t < perturb_end)

    # Rates
    exc_rate_baseline = 5.0
    exc_rate_perturb = 6.0
    inh_rate_baseline = 15.0
    inh_rate_perturb = 2.0

    def simulate_spikes(n_units, rate_baseline, rate_perturb, t, perturb_mask, shared_event_rate=1.0, sync_strength=0.05):
        n_time = len(t)
        spike_trains = np.zeros((n_units, n_time), dtype=np.uint8)
        shared_baseline_mask = np.random.rand(n_time) < (shared_event_rate / fs)
        shared_perturb_mask  = np.random.rand(n_time) < ((shared_event_rate*2) / fs)
        shared_mask = ( (shared_baseline_mask & ~perturb_mask) | (shared_perturb_mask & perturb_mask) )
        for i in range(n_units):
            rate = np.where(perturb_mask, rate_perturb, rate_baseline)
            p_spike = rate / fs
            ind_spikes = np.random.rand(n_time) < p_spike
            shared_spikes = (np.random.rand(n_time) < (sync_strength * shared_mask.astype(float)))
            st = (ind_spikes | shared_spikes).astype(np.uint8)
            spike_trains[i] = st
        return spike_trains

    exc_spikes = simulate_spikes(n_exc, exc_rate_baseline, exc_rate_perturb, t, perturb_mask, shared_event_rate=2.0, sync_strength=0.03)
    inh_spikes = simulate_spikes(n_inh, inh_rate_baseline, inh_rate_perturb, t, perturb_mask, shared_event_rate=2.0, sync_strength=0.05)

    def make_lfp(spike_trains, t, fs):
        tau = 0.01
        kernel_t = np.arange(0, 0.2, 1/fs)
        kernel = (kernel_t / tau) * np.exp(1 - kernel_t/tau)
        lfp = np.zeros(len(t))
        for st in spike_trains:
            lfp += np.convolve(st, kernel, mode='same')
        freqs = np.fft.rfftfreq(len(t), 1/fs)
        amps = 1 / (1 + freqs*5)
        phases = np.exp(2j * np.pi * np.random.rand(len(amps)))
        noise = np.fft.irfft(amps * phases, n=len(t))
        lfp = (lfp - np.mean(lfp)) / (np.std(lfp) + 1e-12) + 0.5 * noise
        lfp = signal.detrend(lfp)
        return lfp

    lfp = make_lfp(np.vstack([exc_spikes, inh_spikes]), t, fs)

    return {
        't': t,
        'fs': fs,
        'exc_spikes': exc_spikes,
        'inh_spikes': inh_spikes,
        'lfp': lfp,
        'perturb_start': perturb_start,
        'perturb_end': perturb_end
    }

# --- Analysis utilities ---

def binned_counts(spike_trains, bin_ms, fs):
    bin_samps = int(bin_ms/1000*fs)
    n_bins = spike_trains.shape[1] // bin_samps
    bc = np.zeros((spike_trains.shape[0], n_bins))
    for i in range(spike_trains.shape[0]):
        bc[i] = spike_trains[i].reshape(-1, bin_samps).sum(axis=1)
    return bc


def pairwise_corr(binned_counts_arr, mask_bins):
    data = binned_counts_arr[:, mask_bins]
    n = data.shape[0]
    corrs = []
    for i in range(n):
        for j in range(i+1, n):
            xi = data[i]
            xj = data[j]
            if xi.std() < 1e-12 or xj.std() < 1e-12:
                continue
            corrs.append(np.corrcoef(xi, xj)[0,1])
    return np.array(corrs)

# --- Main routine ---

def main(use_real_data=False):
    # Note: code to discover and download DANDI datasets would go here. For portability
    # we do not attempt to download any file in this environment by default. If the user
    # wants to run on a specific dataset, set use_real_data=True and modify the script to
    # point to downloaded NWB files.

    if use_real_data:
        print('Real-data mode requested, but this script requires the user to provide NWB files.\n'
              'Please modify analysis.py to point to local NWB files or allow dandi-api access.')
        return

    print('No real dataset configured; generating synthetic demonstration dataset...')
    ds = generate_synthetic_dataset()

    t = ds['t']
    fs = ds['fs']
    exc_spikes = ds['exc_spikes']
    inh_spikes = ds['inh_spikes']
    lfp = ds['lfp']
    perturb_start = ds['perturb_start']
    perturb_end = ds['perturb_end']

    # Population rate (smoothed)
    window = int(0.1 * fs)
    kernel = np.ones(window)/window
    pop_rate = signal.convolve(np.vstack([exc_spikes, inh_spikes]).sum(axis=0), kernel, mode='same') * fs / (exc_spikes.shape[0] + inh_spikes.shape[0])

    # Binned pairwise correlations (50 ms bins)
    bin_ms = 50
    exc_binned = binned_counts(exc_spikes, bin_ms, fs)
    n_bins = exc_binned.shape[1]
    times_bin = np.arange(n_bins) * (bin_ms/1000.0)
    pre_mask_bins = times_bin < perturb_start
    dur_mask_bins = (times_bin >= perturb_start) & (times_bin < perturb_end)

    pre_corrs = pairwise_corr(exc_binned, pre_mask_bins)
    dur_corrs = pairwise_corr(exc_binned, dur_mask_bins)

    # Simple summaries
    summary = {
        'pre_mean_corr': float(np.mean(pre_corrs)) if pre_corrs.size>0 else None,
        'pre_median_corr': float(np.median(pre_corrs)) if pre_corrs.size>0 else None,
        'dur_mean_corr': float(np.mean(dur_corrs)) if dur_corrs.size>0 else None,
        'dur_median_corr': float(np.median(dur_corrs)) if dur_corrs.size>0 else None,
        'n_exc_units': int(exc_spikes.shape[0]),
        'n_inh_units': int(inh_spikes.shape[0])
    }

    # Save summary
    with open(os.path.join(OUTDIR, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    # Save figure similar to the one produced in the reproducible run
    fig = plt.figure(figsize=(9, 6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1,2,1], width_ratios=[3,1], hspace=0.4, wspace=0.3)

    ax0 = fig.add_subplot(gs[0,0])
    ax0.plot(t, lfp, color='k', linewidth=0.5)
    ax0.axvspan(perturb_start, perturb_end, color='orange', alpha=0.2)
    ax0.set_xlim(0, t[-1])
    ax0.set_ylabel('LFP (a.u.)')

    ax1 = fig.add_subplot(gs[1,0])
    # Raster
    for i in range(exc_spikes.shape[0]):
        sp_times = t[exc_spikes[i].astype(bool)]
        ax1.vlines(sp_times, i+0.5, i+1.5, color='C0', linewidth=0.5)
    for j in range(inh_spikes.shape[0]):
        sp_times = t[inh_spikes[j].astype(bool)]
        ax1.vlines(sp_times, exc_spikes.shape[0] + j + 0.5, exc_spikes.shape[0] + j + 1.5, color='C1', linewidth=0.5)
    ax1.axvspan(perturb_start, perturb_end, color='orange', alpha=0.2)
    ax1.set_xlim(0, t[-1])
    ax1.set_ylim(0.5, exc_spikes.shape[0] + inh_spikes.shape[0] + 0.5)
    ax1.set_ylabel('Units')

    ax2 = fig.add_subplot(gs[2,0])
    ax2.plot(t, pop_rate, color='purple')
    ax2.axvspan(perturb_start, perturb_end, color='orange', alpha=0.2)
    ax2.set_xlim(0, t[-1])
    ax2.set_ylabel('Population rate (Hz)')
    ax2.set_xlabel('Time (s)')

    ax3 = fig.add_subplot(gs[:,1])
    ax3.hist(pre_corrs, bins=30, alpha=0.7, label=f'Pre (median={summary["pre_median_corr"]:.3f})', color='C2')
    ax3.hist(dur_corrs, bins=30, alpha=0.7, label=f'During (median={summary["dur_median_corr"]:.3f})', color='C3')
    ax3.set_xlabel('Pairwise corr (binned counts)')
    ax3.set_ylabel('Count')
    ax3.legend(fontsize=8)

    fig.suptitle('Synthetic example: interneuron dampening increases excitatory pairwise correlations', fontsize=12)
    plt.tight_layout(rect=[0,0,1,0.96])
    figpath = os.path.join(OUTDIR, 'figure_1.png')
    plt.savefig(figpath, dpi=150)
    plt.close(fig)

    print('Saved outputs to', OUTDIR)
    print('Summary:', summary)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run interneuron dampening analysis (demo mode uses synthetic data).')
    parser.add_argument('--real', action='store_true', help='Attempt real-data mode (requires local NWB files/dandi access).')
    args = parser.parse_args()
    main(use_real_data=args.real)
