# Depth-Aware Receiver Pre-Configuration for Underwater Acoustic Sensor Networks

Time-varying deep-sea acoustic channel estimation and **depth-aware receiver pre-configuration** for fixed underwater acoustic sensor network (UASN) nodes deployed across thermocline layers.

## Scope Change (v2.0)

> **From**: Depth-Adaptive Variable-Parameter Receiver for AUVs  
> **To**: Depth-Aware Receiver Pre-Configuration for Underwater Acoustic Sensor Networks (UASN)

Four fixed-depth sensor nodes deployed across deep-sea thermocline layers. Each node faces a quasi-static multipath channel characterized by its deployment depth. We propose automatic parameter derivation rules that map channel statistics to optimal DFE equalizer settings — enabling offline pre-configuration without continuous channel estimation.

## Project Structure

```
├── window_*fenduan.m              # Sliding-window CIR estimation (raw WAV → h_matrix)
├── Channel_review_*.m             # BER simulation: RLS-DFE vs LMS-DFE, auto-parameter selection
├── TimeVeryCharater.m             # Channel statistics: coherence time/bandwidth, Doppler, PDP
├── auto_parameter_selector.m      # ★ NEW: Automatic DFE parameter derivation engine
├── utils/extract_multipath_by_threshold.m  # ★ NEW: Automatic multipath extraction
├── uniform_vs_depthaware.m        # ★ NEW: Uniform vs depth-aware receiver comparison
├── complexity_analysis.m          # ★ NEW: Computational complexity quantification
├── network_ber_comparison.m       # ★ NEW: Network-level BER summary
├── parameter_config_table.m       # ★ NEW: Deployment config export (JSON/CSV/LaTeX)
├── plot_*.m                       # Visualization utilities
├── all_PDF.m / PDFfangzhen.m      # Statistical distribution analysis
├── meanCIR.m                      # Mean CIR 3D comparison
└── Depth-Adaptive_.../            # Paper LaTeX source and figures
```

## Network Nodes (Measurement Campaign)

| Node | Depth   | Thermocline Layer | Dataset |
|------|---------|-------------------|---------|
| A    | 475.8 m | Upper             | 时变冲激响应475_0.1_3.mat |
| B    | 624.8 m | Within            | 时变冲激响应624_0.1_3.mat |
| C    | 824.8 m | **Lower ★ bottleneck** | 时变冲激响应824_0.1_3.mat |
| D    | 924.8 m | Below             | 时变冲激响应924_0.1_3.mat |

## Pipeline

### 1. Channel Estimation (`window_*fenduan.m`)
- Raw WAV: transmitted probe, received signal, calibration signal
- Downsample to GCD rate → bandpass filter (600–1000 Hz)
- Cross-correlation synchronization → sliding-window FFT-based CIR estimation
- Output: `h_matrix` (74 time slices × 95999 delay taps, Δt = 0.1 s)

### 2. Channel Characterization (`TimeVeryCharater.m`)
- Time autocorrelation → coherence time (e⁻¹ criterion)
- Frequency autocorrelation → coherence bandwidth
- Doppler power spectrum → Doppler spread
- Power delay profile → RMS delay spread
- **Reusable function**: `computeChannelStats(h_matrix, step_duration, fs)`

### 3. Automatic Parameter Selection (`auto_parameter_selector.m`) ★ NEW
- `tau_rms` → derives `eq_fb_taps` (covers RMS delay spread)
- `T_c` (coherence time) → derives `rls_forget` (tracking speed)
- PDP energy distribution → derives `eq_ff_taps` and `ReferenceTap`
- Output: recommended DFE parameters + FLOPs estimates
- **Replaces** the previous manual look-up table approach

### 4. Communication Performance (`Channel_review_*.m`)
- QPSK modulation (4000 sps, fs = 16 kHz)
- Time-varying convolution (block-wise overlap-add)
- RLS-DFE vs LMS-DFE comparison with automatic parameter selection
- Complexity measurement: `tic/toc` + FLOPs recording

## Key Parameters

| Parameter              | Value               |
|------------------------|---------------------|
| Sampling rate          | 16 kHz              |
| Bandwidth              | 600–1000 Hz         |
| Symbol rate (QPSK)     | 4000 sps            |
| CIR time resolution    | 0.1 s               |
| CIR delay span         | 0–500 ms            |
| DFE forward taps       | Auto-derived (30–120) |
| DFE feedback taps      | Auto-derived (30–800) |
| RLS forgetting factor  | Auto-derived (0.97–0.999) |

## Requirements

- MATLAB R2020a or later
- Communications Toolbox
- Signal Processing Toolbox

## Quick Start

```matlab
% 1. Run channel estimation (one per depth)
window_475fenduan

% 2. Run BER simulation with auto-parameter selection
Channel_review_475

% 3. Compare uniform vs depth-aware deployment
uniform_vs_depthaware

% 4. Generate deployment configuration table
parameter_config_table

% 5. View network-level BER summary
network_ber_comparison
```

## Key Changes from v1.0 (AUV Version)

| v1.0 (Rejected) | v2.0 (This version) |
|-----------------|---------------------|
| Manual selection of 2–4 taps per depth | Automatic -20dB multipath extraction |
| Hardcoded eq_fb_taps per depth | Automatic derivation from channel statistics |
| RLS-DFE only | RLS-DFE vs LMS-DFE comparison |
| No complexity analysis | FLOPs + runtime quantification |
| AUV mobility framing | Fixed UASN node framing |
| Look-up table as core method | Look-up table as engineering implementation of derivation rules |

## License

For academic research purposes.
