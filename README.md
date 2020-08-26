# Pitch Prefect

Real-time (~10ms latency) pitch detection of monophonic signals.

Partly based on the MPM peak picking algorithm from "A Smarter Way to Find
Pitch" by McLeod/Wyvill, but it also analyses the signal spectrogram to further
refine the result. The general problem with MPM is that its pitch estimate is
sometimes off by a factor of `N:1` or `1:N` where `N` is an integer; so it
seems prone to pick up on harmonics, and sometimes even sub-harmonics. MPM has
some fine-tuning parameters, but none of them really solves the problem.
However, harmonics in a monophonic signal are equally spaced in frequency, so
if the fundamental frequency is 100hz, you'll also see spectrogram peaks at
200hz, 300hz, 400hz, and so on. This is used to guess the most likely `N:1` or
`1:N` factor that the MPM pitch estimate is off by, and correct it. The
spectrogram is not directly used to estimate pitch, only the factor. This is
because the spectrogram has a rather low resolution (in order to be real-time).

Placed into the public domain / cc0 licensed. Uses old FFTPACK Fortran code
ported to C, and placed into public domain by Monty of xiph.org (see
http://www.netlib.org/fftpack/fft.c)

The sub directories contain various demos and tools using Pitch Prefect.
