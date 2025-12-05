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
use pxfm::{f_cospi, f_cospif};

mod blackman;
mod hamming;
mod hann;
mod mla;

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
}

pub(crate) trait Trigonometry<V> {
    fn cospi(self) -> V;
}

impl Trigonometry<f32> for f32 {
    #[inline(always)]
    fn cospi(self) -> f32 {
        f_cospif(self)
    }
}

impl Trigonometry<f64> for f64 {
    #[inline(always)]
    fn cospi(self) -> f64 {
        f_cospi(self)
    }
}
