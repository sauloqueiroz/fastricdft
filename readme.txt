--------------- DESCRIPTION
This is a case study of the RIC DFT algorithm proposed in [1] consisting in performing the spectral analysis of the A440 key of a piano.
As widely known, that key concentrates energy on harmonic frequencies multiples of the fundamental frequency 440 Hz, namely, 440 Hz, 880 Hz, 1320 Hz, and so on. Moreover, it is known that the amplitude is inversely proportional to the value of the Harmonic frequency, as one can see in N-point FFT curve given in red across the different plotting files. Our point is that one can "see" the relevant harmonic frequencies with a compressed version of the original signal, as achieved by our C-point RIC DFT algorithm (C in [2, N/2]). As explained in [1], the compression from N to C points is achieved at the expense of only N-1 complex additions (no complex multiplication). Also, we show that a DFT on the resulting C-point signal is equivalent to C coefficients of the full N-point DFT. For example, for C=N/2, the DFT of the compressed signal gives the even-indexed coefficients of the N-point DFT. Because of this equivalence, the C-point compressed signal can be computed by any DFT implementation with no need to algorithmic modifications. Thus, by relying on the C-point signal, one can achieve a spectrum similar to the N-point FFT at a lower computational complexity. In our case study, a very close representation of the spectrum was achieved with the C=N/4-point signal with a complexity of O(ClogC). This demonstrates that our RIC DFT algorithm can capture the trade-off between complexity and spectrum resolution, enabling true improvement on the complexity of the original N-point FFT without sacrificing relevant spectral information.

[1] https://arxiv.org/abs/2504.12551.

--------------- TECHNICAL INFORMATION 
Each file correspond to a different value of C.

Open Octave,
Move to the directory where the files are,
execute source "fastharmonicrics.m"

This will produce different curves of spectral analysis for the audio sample piano_A4_Fs44100Hz.wav.

To get the compacted C-point signal xhat from the original full N-point signal x (N is a power of two), do:
xhat=compactric(x,C); % C is in [2, N/2].

* If this work is useful to you, please cite:
 ********* S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete
           fourier transform square index coefficients,” IEEE Signal Process. Mag.
           (Tips & Tricks), 2025, accepted for publication.
********** Saulo Queiroz, João P. Vilela, Benjamin Koon Kei Ng, Chan-Tong Lam, Edmundo Monteiro., “Fast computation of the discrete
           fourier transform rectangular index coefficients,”, under review by IEEE Signal Processing Letters, avalailable online in https://arxiv.org/abs/2504.12551, 2025.
 *
