/*
 * // Copyright (c) Radzivon Bartoshyk 12/2025. All rights reserved.
 * //
 * // Redistribution and use in source and binary forms, with or without modification,
 * // are permitted provided that the following conditions are met:
 * //
 * // 1.  Redistributions of source code must retain the above copyright notice, this
 * // list of conditions and the following disclaimer.
 * //
 * // 2.  Redistributions in binary form must reproduce the above copyright notice,
 * // this list of conditions and the following disclaimer in the documentation
 * // and/or other materials provided with the distribution.
 * //
 * // 3.  Neither the name of the copyright holder nor the names of its
 * // contributors may be used to endorse or promote products derived from
 * // this software without specific prior written permission.
 * //
 * // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * // DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * // OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
use crate::blackman::blackman_impl;
use crate::hamming::hamming_impl;
use crate::hann::hann_impl;
use crate::kaiser::kaiser_impl;
use crate::slepian::slepian_window;
use num_traits::{AsPrimitive, Float, MulAdd, Signed};
use pxfm::{f_cospi, f_cospif, f_i0, f_i0f};
use std::ops::{AddAssign, Div, Mul, MulAssign, Sub};

mod blackman;
mod hamming;
mod hann;
mod kaiser;
mod mla;
mod slepian;

trait WindowSample:
    Copy
    + Mul<Output = Self>
    + Div<Output = Self>
    + Signed
    + Sub<Output = Self>
    + Float
    + 'static
    + Trigonometry<Self>
    + MulAdd<Self, Output = Self>
    + AddAssign<Self>
    + Float
    + MulAssign<Self>
{
}

impl WindowSample for f64 {}
impl WindowSample for f32 {}

/// Periodic variant of a symmetric window: the symmetric window of length
/// `len + 1` truncated to `len` samples, exactly how SciPy builds `sym=False`
/// windows. A single sample is `[1]`.
fn periodic<V: WindowSample>(len: usize, symmetric: impl Fn(usize) -> Vec<V>) -> Vec<V>
where
    f64: AsPrimitive<V>,
{
    assert!(len > 0, "Windows of size 0 is not defined");
    if len == 1 {
        return vec![1f64.as_()];
    }
    let mut w = symmetric(len + 1);
    w.truncate(len);
    w
}

/// Pxwindow provides methods to generate common window functions used in signal processing.
pub struct Pxwindow {}

impl Pxwindow {
    /// Generates a Hann window of length len in f32 precision.
    pub fn hann_f32(len: usize) -> Vec<f32> {
        hann_impl(len)
    }

    /// Generates a Hann window of length len in f64 precision.
    pub fn hann_f64(len: usize) -> Vec<f64> {
        hann_impl(len)
    }

    /// Generates a Hamming window of length `len` in `f32` precision.
    pub fn hamming_f32(len: usize) -> Vec<f32> {
        hamming_impl(len)
    }

    /// Generates a Hamming window of length `len` in `f64` precision.
    pub fn hamming_f64(len: usize) -> Vec<f64> {
        hamming_impl(len)
    }

    /// Generates a Blackman window of length `len` in `f32` precision.
    pub fn blackman_f32(len: usize) -> Vec<f32> {
        blackman_impl(len)
    }

    /// Generates a Blackman window of length `len` in `f64` precision.
    pub fn blackman_f64(len: usize) -> Vec<f64> {
        blackman_impl(len)
    }

    /// Generates a periodic ("DFT-even") Hann window of length `len` in `f32` precision.
    ///
    /// Equivalent to `scipy.signal.get_window('hann', len)` / `hann(len, sym=False)`:
    /// the symmetric window of length `len + 1` without its last sample. Use it for
    /// spectral analysis (Welch, STFT); the symmetric window is what filter design wants.
    pub fn hann_periodic_f32(len: usize) -> Vec<f32> {
        periodic(len, hann_impl)
    }

    /// Generates a periodic ("DFT-even") Hann window of length `len` in `f64` precision.
    /// See [`Pxwindow::hann_periodic_f32`].
    pub fn hann_periodic_f64(len: usize) -> Vec<f64> {
        periodic(len, hann_impl)
    }

    /// Generates a periodic Hamming window of length `len` in `f32` precision
    /// (`scipy.signal.get_window('hamming', len)`).
    pub fn hamming_periodic_f32(len: usize) -> Vec<f32> {
        periodic(len, hamming_impl)
    }

    /// Generates a periodic Hamming window of length `len` in `f64` precision.
    pub fn hamming_periodic_f64(len: usize) -> Vec<f64> {
        periodic(len, hamming_impl)
    }

    /// Generates a periodic Blackman window of length `len` in `f32` precision
    /// (`scipy.signal.get_window('blackman', len)`).
    pub fn blackman_periodic_f32(len: usize) -> Vec<f32> {
        periodic(len, blackman_impl)
    }

    /// Generates a periodic Blackman window of length `len` in `f64` precision.
    pub fn blackman_periodic_f64(len: usize) -> Vec<f64> {
        periodic(len, blackman_impl)
    }

    /// Generates a Slepian window of length `len` in `f32` precision.
    pub fn slepian_f32(len: usize, nw: f64) -> Vec<f32> {
        slepian_window(len, nw).iter().map(|x| *x as f32).collect()
    }

    /// Generates a Slepian window of length `len` in `f64` precision.
    pub fn slepian_f64(len: usize, nw: f64) -> Vec<f64> {
        slepian_window(len, nw)
    }

    /// Generates a Kaiser window of length `len` in `f32` precision.
    pub fn kaiser_f32(len: usize, beta: f32) -> Vec<f32> {
        kaiser_impl(len, beta)
    }

    /// Generates a Slepian window of length `len` in `f64` precision.
    pub fn kaiser_f64(len: usize, beta: f64) -> Vec<f64> {
        kaiser_impl(len, beta)
    }
}

pub(crate) trait Trigonometry<V> {
    fn cospi(self) -> V;
    fn i0(self) -> V;
}

impl Trigonometry<f32> for f32 {
    #[inline(always)]
    fn cospi(self) -> f32 {
        f_cospif(self)
    }

    #[inline(always)]
    fn i0(self) -> f32 {
        f_i0f(self)
    }
}

impl Trigonometry<f64> for f64 {
    #[inline(always)]
    fn cospi(self) -> f64 {
        f_cospi(self)
    }

    #[inline(always)]
    fn i0(self) -> f64 {
        f_i0(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(got: &[f64], expected: &[f64]) {
        assert_eq!(got.len(), expected.len());
        for (g, e) in got.iter().zip(expected) {
            assert!((g - e).abs() < 1e-12, "got {g}, expected {e}");
        }
    }

    /// Reference: `scipy.signal.windows.{hann,hamming,blackman}(n, sym=False)`.
    #[test]
    fn periodic_windows_match_scipy() {
        assert_close(
            &Pxwindow::hann_periodic_f64(8),
            &[
                0.0,
                0.14644660940672627,
                0.5,
                0.8535533905932737,
                1.0,
                0.8535533905932737,
                0.5,
                0.14644660940672627,
            ],
        );
        assert_close(
            &Pxwindow::hann_periodic_f64(7),
            &[
                0.0,
                0.18825509907063326,
                0.6112604669781572,
                0.9504844339512095,
                0.9504844339512095,
                0.6112604669781572,
                0.18825509907063326,
            ],
        );
        assert_close(
            &Pxwindow::hamming_periodic_f64(8),
            &[
                0.08000000000000007,
                0.21473088065418822,
                0.54,
                0.865269119345812,
                1.0,
                0.865269119345812,
                0.54,
                0.21473088065418822,
            ],
        );
        assert_close(
            &Pxwindow::blackman_periodic_f64(7),
            &[
                -1.3877787807814457e-17,
                0.09045342435412808,
                0.45918295754596367,
                0.9203636180999082,
                0.9203636180999082,
                0.45918295754596367,
                0.09045342435412808,
            ],
        );
        assert_eq!(Pxwindow::hann_periodic_f64(1), vec![1.0]);
        assert_eq!(Pxwindow::hann_periodic_f32(4), vec![0.0, 0.5, 1.0, 0.5]);
    }
}
