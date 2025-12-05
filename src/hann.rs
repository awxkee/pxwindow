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
use crate::Trigonometry;
use num_traits::{AsPrimitive, Float, Signed};
use std::ops::{Div, Mul, Sub};

#[inline(always)]
pub(crate) fn hann_impl<
    V: Copy
        + Mul<Output = V>
        + Div<Output = V>
        + Signed
        + Sub<Output = V>
        + Float
        + 'static
        + Trigonometry<V>,
>(
    len: usize,
) -> Vec<V>
where
    f64: AsPrimitive<V>,
    usize: AsPrimitive<V>,
{
    assert!(len > 0, "Windows of size 0 is not defined");

    let mut w = vec![V::zero(); len];

    if len == 1 {
        return vec![1f64.as_()];
    }

    let denom: V = (1. / (len - 1) as f64).as_();

    for (n, dst) in w.iter_mut().enumerate() {
        let r: V = (2.0f64.as_() * n.as_()) * denom;
        let v = 0.5f64.as_() * (1.0f64.as_() - r.cospi());
        *dst = v;
    }

    w
}
