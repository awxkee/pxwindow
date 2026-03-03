/*
 * // Copyright (c) Radzivon Bartoshyk 3/2026. All rights reserved.
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

pub(crate) fn slepian_window(n: usize, nw: f64) -> Vec<f64> {
    let w = nw / n as f64;
    let mut diag = vec![0.0f64; n];
    let mut off = vec![0.0f64; n - 1];
    for (i, dst) in diag.iter_mut().enumerate() {
        *dst = ((n as f64 - 1.0) / 2.0 - i as f64).powi(2) * (2.0 * w).cospi();
    }
    for (i, dst) in off.iter_mut().enumerate().take(n - 1) {
        let k = (i + 1) as f64;
        *dst = k * (n as f64 - k) * 0.5;
    }
    let lambda = largest_eigenvalue(&diag, &off, n);
    let mut v = inverse_iteration(&diag, &off, n, lambda);
    if v.iter().sum::<f64>() < 0.0 {
        for x in &mut v {
            *x = -*x;
        }
    }
    let max_val = 1. / v.iter().cloned().fold(0.0f64, |acc, x| acc.max(x.abs()));
    for x in &mut v {
        *x *= max_val;
    }
    if n.is_multiple_of(2) {
        let m2 = (n as f64).powi(2);
        let correction = m2 / (m2 + nw);
        for x in &mut v {
            *x *= correction;
        }
    }
    v
}

fn sturm_count_above(diag: &[f64], off: &[f64], n: usize, mu: f64) -> usize {
    let mut below = 0usize;
    let mut d = diag[0] - mu;
    if d < 0.0 {
        below += 1;
    }
    for i in 1..n {
        let denom = if d.abs() < 1e-300 {
            1e-300_f64.copysign(d)
        } else {
            d
        };
        d = (diag[i] - mu) - off[i - 1].powi(2) / denom;
        if d < 0.0 {
            below += 1;
        }
    }
    n - below
}

fn gershgorin(diag: &[f64], off: &[f64], n: usize) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..n {
        let r = if i == 0 {
            off[0].abs()
        } else if i == n - 1 {
            off[n - 2].abs()
        } else {
            off[i - 1].abs() + off[i].abs()
        };
        lo = lo.min(diag[i] - r);
        hi = hi.max(diag[i] + r);
    }
    (lo - 1.0, hi + 1.0)
}

fn largest_eigenvalue(diag: &[f64], off: &[f64], n: usize) -> f64 {
    let (mut lo, mut hi) = gershgorin(diag, off, n);
    for _ in 0..200 {
        let mid = (lo + hi) / 2.0;
        if sturm_count_above(diag, off, n, mid) >= 1 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-14 {
            break;
        }
    }
    (lo + hi) / 2.0
}

fn tridiag_solve(diag: &[f64], off: &[f64], n: usize, shift: f64, b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.0f64; n - 1];
    let mut d = vec![0.0f64; n];
    let mut w = diag[0] - shift;
    if w.abs() < 1e-14 {
        w = 1e-14;
    }
    d[0] = b[0] / w;
    if n > 1 {
        c[0] = off[0] / w;
    }
    for i in 1..n {
        let mut denom = (diag[i] - shift) - off[i - 1] * c[i - 1];
        if denom.abs() < 1e-14 {
            denom = 1e-14;
        }
        d[i] = (b[i] - off[i - 1] * d[i - 1]) / denom;
        if i < n - 1 {
            c[i] = off[i] / denom;
        }
    }
    let mut x = vec![0.0f64; n];
    x[n - 1] = d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d[i] - c[i] * x[i + 1];
    }
    x
}

fn inverse_iteration(diag: &[f64], off: &[f64], n: usize, lambda: f64) -> Vec<f64> {
    let mut v: Vec<f64> = (0..n)
        .map(|i| (i as f64 / (n as f64 - 1.0)).cospi())
        .collect();
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in &mut v {
        *x /= norm;
    }
    let shift = lambda + 1e-8;
    for _ in 0..100 {
        let w2 = tridiag_solve(diag, off, n, shift, &v);
        let nw: f64 = w2.iter().map(|x| x * x).sum::<f64>().sqrt();
        if nw < 1e-300 {
            break;
        }
        let v_new: Vec<f64> = w2.iter().map(|x| x / nw).collect();
        let diff = [1.0_f64, -1.0]
            .iter()
            .map(|&s| {
                v_new
                    .iter()
                    .zip(v.iter())
                    .map(|(a, b)| (a - s * b).abs())
                    .fold(0.0f64, f64::max)
            })
            .fold(f64::INFINITY, f64::min);
        v = v_new;
        if diff < 1e-13 {
            break;
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn length() {
        assert_eq!(slepian_window(64, 4.0 / 64.0).len(), 64);
    }

    #[test]
    fn symmetric() {
        let w = slepian_window(64, 4.0 / 64.0);
        let n = w.len();
        for i in 0..n / 2 {
            assert!(close(w[i], w[n - 1 - i], 1e-12), "not symmetric at {i}");
        }
    }

    #[test]
    fn values_bounded() {
        let w = slepian_window(64, 4.0 / 64.0);
        for &v in &w {
            assert!(v.abs() <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn all_positive_k0() {
        // k=0 DPSS has no sign changes
        let w = slepian_window(64, 4.0 / 64.0);
        assert!(w.iter().all(|&v| v > 0.0));
    }

    #[test]
    fn matches_scipy_nw4_n8() {
        let expected = [
            0.01012579, 0.08131032, 0.29892319, 0.65967608, 0.96246391, 0.96246391, 0.65967608,
            0.29892319, 0.08131032, 0.01012579,
        ];
        let w = slepian_window(10, 3.9);
        for (got, &exp) in w.iter().zip(expected.iter()) {
            assert!(close(*got, exp, 1e-6), "got {got:.7} expected {exp:.7}");
        }
    }
}
