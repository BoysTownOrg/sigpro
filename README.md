SIGPRO

Synopsis

SIGPRO is a library of signal-processing functions designed to assist in the development of auditory research software. Current functions include random number generators, fft, inverse fft, frequency shaping (filtering), and sample rate conversion. Limited support is provided for loading and saving MATLAB (version 4) binary (MAT) files. 

Code Example

    // Butterworth filter
    sp_butter(b, a, no, wn, ft);        // filter coefficients
    sp_filter(b, nc, a, nc, x, y, np);  // impulse response

Motivation

To provide C programs with basic signal-processing functions available in MATLAB.

Installation

Download repo from https://github.com/BTNRH/sigpro. Makefiles are provided for building test programs at Linux, MacOS, or MinGW command lines. A solution file is provided in the VS9 folder for building under Visual Studio.

API Reference

The API is described in the User Manual at https://github.com/BTNRH/sigpro/blob/master/sigpro.pdf.

Tests

Test programs are provided to demonstrate several features of functions included in the library:

    tst_afd - test analog-filter design
    tst_cpx - test writing & reading complex data in MAT files
    tst_fft - test real & complex FFT and inverse FFT
    tst_mat - test MAT file save & load
    tst_min - test fmins function
    tst_shp - test frequency-shaping functions
    tst_src - test sampling-rate conversion
    tst_tb  - test tone-burst waveform calculation
    tst_wav - test WAV file read & write
    tst_xfr - test transfer-function computation

Contributors

Report bugs to Stephen.Neely@boystown.org.

License

Creative Commons?

