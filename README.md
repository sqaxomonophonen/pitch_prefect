# Pitch Prefect

Real-time (~10ms latency) pitch detection of monophonic signals.

It's based on the MPM algorithm from "A Smarter Way to Find Pitch" by McLeod
and Wyvill, but it uses spectrogram analysis to further refine the result.

The MPM algorithm produces accurate `f0` (fundamental frequency / pitch)
estimates most of the time, but is sometimes off by a factor of `N:1` or `1:N`
where `N` is an integer; it means that it's prone to pick up on harmonics, and
sometimes even sub-harmonics (I'm not sure how that makes sense, but it does
happen). MPM has some fine-tuning parameters, but none of them really solves
the problem. MPM uses a NSDF (Normalized Square Distance Function) which is
similar to auto-correlation, but looking at NSDF plots suggests that the
problem can't easily be solved with the NSDF alone; it outputs maxima at delays
where the "self-similarity" is highest, and these delays can be directly
translated into frequencies, but sometimes the best candidate (highest maximum)
corresponds a harmonic or sub-harmonic rather than `f0`.

However, harmonics in a monophonic signal are equally spaced in frequency; if
`f0` is 100hz, you'll also see spectrogram peaks at 200hz, 300hz, 400hz, and so
on. If the MPM error factor is `2:1` or `3:1`, it'll report `f0=200hz` or
`f0=300hz` rather than the correct answer (`f0=100hz`). The spectogram analysis
goes through a static list of likely `1:N` and `N:1` candidates and finds the
best match. The spectogram analysis does not attempt to guess `f0` on its own
because the spectogram is rather low in resolution (a real-time requirement).


# Issues

 - Alpha software: it hasn't seen any heavy duty usage yet

 - `f0` estimates are noisy when signal-to-noise is low. This could probably be
   solved with filtering (at the expense of latency of course)

 - There's room for optimizations! NSDF calculation could probably be
   vectorized (SIMD)


# Examples

The sub directories contain various demos and tools using Pitch Prefect.


# License

Placed into the public domain / cc0 licensed. Uses old FFTPACK Fortran code
ported to C, and placed into public domain by Monty of xiph.org (see
http://www.netlib.org/fftpack/fft.c)
