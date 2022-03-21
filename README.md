# OSC-DECOMP (oscillator decomposition of time series data)

This is a statistical method for extracting oscillators from (univariate/multivariate) time series data [1,2].

The main functions are osc_decomp, osc_plot, osc_phase_plot and osc_spectrum_plot.
Please see demo_uni.m and demo_multi.m for their usage.

- demo_uni.m is the script for decomposing the Canadian Lynx data [1].

- demo_multi.m is the script for decomposing the north/south sunspot data [2].

There are several updates from the original papers [1,2] in implementation (optimization, confidence interval [3]).
Please see detail.pdf for their details.

If you have any comments or bug reports, please contact Takeru Matsuda (takeru.matsuda@riken.jp).

References:

[1] T. Matsuda and F. Komaki. Time series decomposition into oscillation components and phase estimation. Neural Computation, Vol. 29, pp. 332--367, 2017. 

[2] T. Matsuda and F. Komaki. Multivariate time series decomposition into oscillation components. Neural Computation, Vol. 29, pp. 2055--2075, 2017. 

[3] T. Matsuda, F. Homae, H. Watanabe, G. Taga and F. Komaki. Oscillator decomposition of infant fNIRS data. PLOS Computational Biology, 2022.
