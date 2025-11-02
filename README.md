# Rayleigh–Rician Channel Analysis

This project presents a **simulation and analysis of Rayleigh and Rician fading channels** and investigates their **BER performance** under varying **SNR** and **time variation** levels.

## Files

- **rayleigh_rician_fading_analysis_report.pdf**  
  Final report describing the theoretical background, channel modeling, simulation setup, and results.
- **rayleigh_rician_fading_analysis_code.m**  
  MATLAB code that simulates multiple channel types (Rayleigh, Rician, AWGN) and plots impulse responses and BER comparisons.

## Overview

The simulation covers:
- Flat and frequency-selective fading models  
- Time-invariant and time-varying cases  
- Effects of different SNR and Doppler shift values  
- Comparative BER plots across channel types  

## Results

The analysis demonstrates that:
- Rician channels generally outperform Rayleigh in moderate SNRs due to LOS components.  
- Increasing Doppler shift (Tm) increases BER due to faster time variation.  
- Flat fading channels are less dispersive but suffer deeper fades.  

## Usage

Run the MATLAB file to reproduce the results:
```matlab
rayleigh_rician_fading_analysis_code
