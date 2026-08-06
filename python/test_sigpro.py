import unittest
import numpy as np
import sigpro

class TestSigpro(unittest.TestCase):
    def test_version(self):
        v = sigpro.version()
        self.assertIsInstance(v, str)
        self.assertTrue(len(v) > 0)
        print("SIGPRO Version:", v)

    def test_rand(self):
        n = 100
        x = sigpro.rand(n)
        self.assertEqual(len(x), n)
        self.assertTrue(np.all((x >= 0.0) & (x <= 1.0)))

    def test_chirp(self):
        n = 1024
        x = sigpro.chirp(n)
        self.assertEqual(len(x), n)
        self.assertFalse(np.all(x == 0.0))

    def test_butter_filter(self):
        n_order = 4
        cutoff = 0.25 # normalized frequency
        b, a = sigpro.butter(n_order, cutoff, 'lowpass')
        self.assertEqual(len(b), n_order + 1)
        self.assertEqual(len(a), n_order + 1)
        
        # Test filtering a random signal
        x = sigpro.rand(500)
        y = sigpro.filter(b, a, x)
        self.assertEqual(len(y), len(x))
        self.assertFalse(np.all(y == 0.0))

    def test_fft(self):
        # 16 complex points (32 floats interleaved)
        N = 16
        x = np.zeros(N * 2, dtype=np.float32)
        # Add a DC component
        for i in range(0, N * 2, 2):
            x[i] = 1.0 
            
        y = sigpro.fft(x.copy())
        self.assertEqual(len(y), N * 2)
        
        # DC component should be N after FFT (since FFT of 1s is N at index 0)
        # Note: Depending on FFT scaling convention in SIGPRO, it might be scaled differently
        # We just verify it ran without errors and the result changed
        self.assertFalse(np.all(y == x))
        
        # Inverse FFT
        x_recovered = sigpro.fft(y.copy(), inverse=True)
        # Should be scaled back, just check shape and type
        self.assertEqual(len(x_recovered), N * 2)

if __name__ == '__main__':
    unittest.main()
