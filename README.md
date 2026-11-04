# OSC-DECOMP: MATLAB Code for Oscillator Decomposition

A MATLAB implementation of **Oscillator Decomposition (OSC-DECOMP)** for time series analysis.

## Requirements

- **MATLAB**
- **Optimization Toolbox**

## Usage

Check the following demo scripts to see how to use the functions:
- `demo_lynx.m`
- `demo_sunspot.m`
- `demo_eeg.m`

### Main Functions
- `osc_decomp`
- `osc_plot`
- `osc_phase_plot`

> **Note:** There are several updates from the original papers [[1](#references), [2](#references)] in implementation (e.g., optimization algorithms and confidence intervals [[3](#references)]).

## References

1. T. Matsuda and F. Komaki. Time series decomposition into oscillation components and phase estimation. *Neural Computation*, Vol. 29, pp. 332–367, 2017.
2. T. Matsuda and F. Komaki. Multivariate time series decomposition into oscillation components. *Neural Computation*, Vol. 29, pp. 2055–2075, 2017.
3. T. Matsuda, F. Homae, H. Watanabe, G. Taga and F. Komaki. Oscillator decomposition of infant fNIRS data. *PLOS Computational Biology*, 2022.

## Contact

For any comments, questions, or bug reports, please contact Takeru Matsuda ([matsuda@mist.i.u-tokyo.ac.jp](mailto:matsuda@mist.i.u-tokyo.ac.jp)).
