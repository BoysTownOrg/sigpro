import ctypes
import numpy as np
import os

# Try to find the shared library
lib_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "libsigpro.so")
if not os.path.exists(lib_path):
    # Try assuming it's a Mac dylib or in path
    lib_path = "libsigpro.so"

try:
    libsigpro = ctypes.CDLL(lib_path)
except OSError:
    raise RuntimeError(f"Could not load libsigpro.so from {lib_path}. Please build it first using 'make libsigpro.so'")

# -------------------------------------------------------------------------
# sp_version
# -------------------------------------------------------------------------
libsigpro.sp_version.restype = ctypes.c_char_p

def version():
    """Returns the SIGPRO library version."""
    return libsigpro.sp_version().decode('utf-8')

# -------------------------------------------------------------------------
# sp_butter
# -------------------------------------------------------------------------
libsigpro.sp_butter.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int
]

def butter(n, wn, btype='lowpass'):
    """
    Butterworth filter design.
    Returns (b, a) coefficients.
    btype can be 'lowpass', 'highpass', 'bandpass', 'bandstop'.
    """
    ft_map = {'lowpass': 0, 'highpass': 1, 'bandpass': 2, 'bandstop': 3}
    ft = ft_map.get(btype.lower(), 0)
    
    wn_arr = np.atleast_1d(wn).astype(np.float32)
    
    b = np.zeros(n + 1, dtype=np.float32)
    a = np.zeros(n + 1, dtype=np.float32)
    
    libsigpro.sp_butter(b, a, n, wn_arr, ft)
    return b, a

# -------------------------------------------------------------------------
# sp_fft / sp_ifft
# -------------------------------------------------------------------------
libsigpro.sp_fft.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int
]
libsigpro.sp_ifft.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int
]

def fft(x, inverse=False):
    """
    Complex-to-complex FFT or iFFT in place.
    x must be a 1D complex array (interleaved real, imag) of size 2*N where N is power of 2.
    Returns the transformed array (which is modified in-place).
    """
    x_c = np.ascontiguousarray(x, dtype=np.float32)
    n = len(x_c) // 2
    if inverse:
        err = libsigpro.sp_ifft(x_c, n)
    else:
        err = libsigpro.sp_fft(x_c, n)
    if err != 0:
        raise ValueError(f"FFT returned error code {err}")
    return x_c

# -------------------------------------------------------------------------
# sp_chirp
# -------------------------------------------------------------------------
libsigpro.sp_chirp.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int
]

def chirp(n):
    """Generates a frequency sweep signal of length n."""
    x = np.zeros(n, dtype=np.float32)
    libsigpro.sp_chirp(x, n)
    return x

# -------------------------------------------------------------------------
# sp_filter
# -------------------------------------------------------------------------
libsigpro.sp_filter.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int
]

def filter(b, a, x):
    """
    IIR filter.
    Returns the filtered signal y.
    """
    b_c = np.ascontiguousarray(b, dtype=np.float32)
    a_c = np.ascontiguousarray(a, dtype=np.float32)
    x_c = np.ascontiguousarray(x, dtype=np.float32)
    y_c = np.zeros_like(x_c)
    
    err = libsigpro.sp_filter(b_c, len(b_c), a_c, len(a_c), x_c, y_c, len(x_c))
    if err != 0:
        raise ValueError(f"sp_filter returned error code {err}")
    return y_c

# -------------------------------------------------------------------------
# sp_rand
# -------------------------------------------------------------------------
libsigpro.sp_rand.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int
]

def rand(n):
    """Generates uniform random values of length n."""
    x = np.zeros(n, dtype=np.float32)
    libsigpro.sp_rand(x, n)
    return x
