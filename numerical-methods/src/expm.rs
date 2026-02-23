//! Matrix exponential computation for CTMC generator matrices.
//!
//! Two interfaces:
//!   1. `expm(a)` — full dense matrix exponential via Padé(13) + scaling-and-squaring.
//!      Use when NM < ~1000 or when you need the full matrix A = exp(G·Δt) for reuse.
//!
//!   2. `expm_krylov_action(matvec, v, n, krylov_dim, tol)` — computes exp(A)·v without
//!      forming exp(A). Use when NM > ~1000 and you only need matrix-vector products.
//!
//! Reference for Padé method:
//!   Higham, N. J. (2005). "The Scaling and Squaring Method for the Matrix Exponential
//!   Revisited." SIAM J. Matrix Anal. Appl. 26(4), 1179-1193.
//!
//! Reference for Krylov:
//!   Saad, Y. (1992). "Analysis of some Krylov subspace approximations to the matrix
//!   exponential operator." SIAM J. Numer. Anal. 29(1), 209-228.

use nalgebra::{DMatrix, DVector};

/// Padé(13) coefficients b_0, ..., b_13 from Higham (2005) Table 1.
const PADE13_COEFFS: [f64; 14] = [
    64764752532480000.0,
    32382376266240000.0,
    7771770303897600.0,
    1187353796428800.0,
    129060195264000.0,
    10559470521600.0,
    670442572800.0,
    33522128640.0,
    1323241920.0,
    40840800.0,
    960960.0,
    16380.0,
    182.0,
    1.0,
];

/// Dense matrix exponential exp(A) via Padé(13) + scaling-and-squaring.
///
/// Algorithm matches MATLAB expm and scipy linalg.expm. Cost is O(n^3)
/// with ~6-10 matrix multiplications.
///
/// # Panics
/// Panics if the input is not square.
pub fn expm(a: &DMatrix<f64>) -> DMatrix<f64> {
    let n = a.nrows();
    assert_eq!(n, a.ncols(), "expm: matrix must be square");

    if n == 0 {
        return DMatrix::zeros(0, 0);
    }
    if n == 1 {
        return DMatrix::from_element(1, 1, a[(0, 0)].exp());
    }

    // Scaling: choose s so that ||A/2^s||_1 <= theta_13.
    let theta_13 = 5.371920351148152_f64;
    let norm_a = matrix_1_norm(a);
    let s = if norm_a <= theta_13 {
        0_u32
    } else {
        ((norm_a / theta_13).log2().ceil() as u32).max(1)
    };

    let a_scaled = a * 2.0_f64.powi(-(s as i32));
    let result = pade13_approximant(&a_scaled);
    square_repeatedly(result, s)
}

/// Compute exp((G - r·I) · dt) — the transition matrix for one timestep.
pub fn expm_generator(g: &DMatrix<f64>, r: f64, dt: f64) -> DMatrix<f64> {
    let n = g.nrows();
    let shifted = (g - DMatrix::identity(n, n) * r) * dt;
    expm(&shifted)
}

fn matrix_1_norm(a: &DMatrix<f64>) -> f64 {
    (0..a.ncols())
        .map(|j| (0..a.nrows()).map(|i| a[(i, j)].abs()).sum::<f64>())
        .fold(0.0_f64, f64::max)
}

/// Padé(13) rational approximant using Higham's Algorithm 2.3.
///
/// Computes the numerator polynomial U and denominator V, then solves
/// (V - U) · R = (V + U) for R ≈ exp(A).
fn pade13_approximant(a: &DMatrix<f64>) -> DMatrix<f64> {
    let n = a.nrows();
    let id = DMatrix::identity(n, n);
    let b = &PADE13_COEFFS;

    let a2 = a * a;
    let a4 = &a2 * &a2;
    let a6 = &a4 * &a2;

    // U = A · (A6·(b13·A6 + b11·A4 + b9·A2) + b7·A6 + b5·A4 + b3·A2 + b1·I)
    let u_high = &a6 * b[13] + &a4 * b[11] + &a2 * b[9];
    let u_poly = &a6 * &u_high + &a6 * b[7] + &a4 * b[5] + &a2 * b[3] + &id * b[1];
    let u = a * &u_poly;

    // V = A6·(b12·A6 + b10·A4 + b8·A2) + b6·A6 + b4·A4 + b2·A2 + b0·I
    let v_high = &a6 * b[12] + &a4 * b[10] + &a2 * b[8];
    let v = &a6 * &v_high + &a6 * b[6] + &a4 * b[4] + &a2 * b[2] + &id * b[0];

    // Solve (V - U) · R = (V + U)
    let lhs = &v - &u;
    let rhs = v + u;
    lhs.lu()
        .solve(&rhs)
        .expect("expm: singular (V-U) in Padé(13) — generator is degenerate")
}

fn square_repeatedly(mut m: DMatrix<f64>, s: u32) -> DMatrix<f64> {
    for _ in 0..s {
        m = &m * &m;
    }
    m
}

/// Compute exp(A) · v without forming exp(A).
///
/// Uses Arnoldi iteration to build the Krylov basis V_k and the
/// upper Hessenberg projection H_k, then evaluates:
///   exp(A)·v ≈ ‖v‖ · V_k · exp(H_k) · e_1
///
/// # Arguments
/// * `matvec` — closure computing A·x for a given x
/// * `v` — input vector of length `n`
/// * `n` — system dimension
/// * `krylov_dim` — maximum Krylov dimension (typically 30-60)
/// * `tol` — breakdown tolerance (typically 1e-10 to 1e-12)
pub fn expm_krylov_action<F>(
    matvec: &F,
    v: &[f64],
    n: usize,
    krylov_dim: usize,
    tol: f64,
) -> Vec<f64>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let beta = vec_norm(v);
    if beta < 1e-30 {
        return vec![0.0; n];
    }

    let k_max = krylov_dim.min(n);
    let mut v_cols: Vec<Vec<f64>> = Vec::with_capacity(k_max + 1);
    // H is stored as (k_max+1) rows × k_max cols (upper Hessenberg)
    let mut h = vec![vec![0.0_f64; k_max]; k_max + 1];

    v_cols.push(v.iter().map(|&x| x / beta).collect());

    let mut k_actual = k_max;
    for j in 0..k_max {
        let mut w = matvec(&v_cols[j]);

        // Modified Gram-Schmidt
        for i in 0..=j {
            let hij = vec_dot(&w, &v_cols[i]);
            h[i][j] = hij;
            for l in 0..n {
                w[l] -= hij * v_cols[i][l];
            }
        }

        let h_next = vec_norm(&w);
        h[j + 1][j] = h_next;

        if h_next < tol * beta {
            k_actual = j + 1;
            break;
        }

        if j + 1 < k_max {
            v_cols.push(w.iter().map(|&x| x / h_next).collect());
        }
    }

    // exp(H_k) · e_1 = first column of exp(H_k); H_k is small, use dense Padé
    let h_small = DMatrix::from_fn(k_actual, k_actual, |i, j| h[i][j]);
    let exp_h = expm(&h_small);

    let mut result = vec![0.0; n];
    for j in 0..k_actual {
        let coeff = beta * exp_h[(j, 0)];
        for l in 0..n {
            result[l] += coeff * v_cols[j][l];
        }
    }
    result
}

/// Krylov action with a dense DMatrix as the operator.
pub fn expm_action_dense(a: &DMatrix<f64>, v: &[f64], krylov_dim: usize) -> Vec<f64> {
    let n = a.nrows();
    let matvec = |x: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|j| a[(i, j)] * x[j]).sum())
            .collect()
    };
    expm_krylov_action(&matvec, v, n, krylov_dim, 1e-12)
}

/// Compressed Sparse Row matrix. G has ~5-6 nonzeros per row so CSR
/// gives O(nnz) mat-vecs vs O(n^2) for dense.
pub struct CsrMatrix {
    pub nrows: usize,
    pub ncols: usize,
    /// row_ptr[i]..row_ptr[i+1] gives the range of entries for row i
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<usize>,
    pub values: Vec<f64>,
}

impl CsrMatrix {
    /// Build CSR from a dense DMatrix, dropping entries below `threshold`.
    pub fn from_dense(a: &DMatrix<f64>, threshold: f64) -> Self {
        let nrows = a.nrows();
        let ncols = a.ncols();
        let mut row_ptr = Vec::with_capacity(nrows + 1);
        let mut col_idx = Vec::new();
        let mut values = Vec::new();

        row_ptr.push(0);
        for i in 0..nrows {
            for j in 0..ncols {
                let v = a[(i, j)];
                if v.abs() > threshold {
                    col_idx.push(j);
                    values.push(v);
                }
            }
            row_ptr.push(col_idx.len());
        }

        Self { nrows, ncols, row_ptr, col_idx, values }
    }

    /// Build CSR from COO triplets. Duplicate (row, col) entries are summed.
    pub fn from_coo(
        nrows: usize,
        ncols: usize,
        rows: &[usize],
        cols: &[usize],
        vals: &[f64],
    ) -> Self {
        // Sort by (row, col)
        let mut entries: Vec<(usize, usize, f64)> = rows
            .iter()
            .zip(cols)
            .zip(vals)
            .map(|((&r, &c), &v)| (r, c, v))
            .collect();
        entries.sort_unstable_by_key(|e| (e.0, e.1));

        // Deduplicate by summing entries at the same (row, col)
        let mut deduped: Vec<(usize, usize, f64)> = Vec::with_capacity(entries.len());
        for (r, c, v) in entries {
            match deduped.last_mut() {
                Some(last) if last.0 == r && last.1 == c => last.2 += v,
                _ => deduped.push((r, c, v)),
            }
        }

        // Build CSR: count entries per row, prefix-sum for row_ptr, fill arrays
        let mut row_ptr = vec![0usize; nrows + 1];
        for &(r, _, _) in &deduped {
            row_ptr[r + 1] += 1;
        }
        for i in 0..nrows {
            row_ptr[i + 1] += row_ptr[i];
        }

        let mut col_idx = Vec::with_capacity(deduped.len());
        let mut values = Vec::with_capacity(deduped.len());
        for (_, c, v) in deduped {
            col_idx.push(c);
            values.push(v);
        }

        Self { nrows, ncols, row_ptr, col_idx, values }
    }

    /// Sparse matrix-vector product y = A · x.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        debug_assert_eq!(x.len(), self.ncols);
        let mut y = vec![0.0; self.nrows];
        for i in 0..self.nrows {
            let mut sum = 0.0;
            for idx in self.row_ptr[i]..self.row_ptr[i + 1] {
                sum += self.values[idx] * x[self.col_idx[idx]];
            }
            y[i] = sum;
        }
        y
    }

    /// Shifted mat-vec y = (A - shift·I) · x without modifying the matrix.
    pub fn matvec_shifted(&self, x: &[f64], shift: f64) -> Vec<f64> {
        let mut y = self.matvec(x);
        for i in 0..self.nrows.min(self.ncols) {
            y[i] -= shift * x[i];
        }
        y
    }

    pub fn nnz(&self) -> usize {
        self.values.len()
    }
}

/// Krylov action for the shifted generator: exp((G - r·I)·dt) · v.
///
/// Use this inside Algorithm 3.1's time loop when NM is large.
pub fn expm_generator_action_sparse(
    g_csr: &CsrMatrix,
    r: f64,
    dt: f64,
    v: &[f64],
    krylov_dim: usize,
) -> Vec<f64> {
    let n = g_csr.nrows;
    let matvec = |x: &[f64]| -> Vec<f64> {
        let mut y = g_csr.matvec_shifted(x, r);
        for val in y.iter_mut() {
            *val *= dt;
        }
        y
    };
    expm_krylov_action(&matvec, v, n, krylov_dim, 1e-12)
}

/// Strategy for computing the matrix exponential action.
pub enum MatExpStrategy {
    /// Precompute A = exp((G-r·I)·dt) once (dense), reuse A·v each step.
    /// Best for NM ≲ 1000. One-time O(n^3), per-step O(n^2).
    DensePade,
    /// Krylov at every step; never form A.
    /// Best for NM ≳ 1000. Per-step O(k·nnz) with k = krylov_dim.
    Krylov { krylov_dim: usize },
}

/// Precomputed operator A for use in Algorithm 3.1's backward loop.
pub enum MatExpOperator {
    Dense { a_matrix: DMatrix<f64> },
    Sparse { g_shifted_dt: CsrMatrix, krylov_dim: usize },
}

impl MatExpOperator {
    /// Build the operator from the combined generator G, risk-free rate r, and timestep dt.
    pub fn new(g: &DMatrix<f64>, r: f64, dt: f64, strategy: MatExpStrategy) -> Self {
        let n = g.nrows();
        match strategy {
            MatExpStrategy::DensePade => {
                MatExpOperator::Dense { a_matrix: expm_generator(g, r, dt) }
            }
            MatExpStrategy::Krylov { krylov_dim } => {
                let shifted = (g - DMatrix::identity(n, n) * r) * dt;
                MatExpOperator::Sparse {
                    g_shifted_dt: CsrMatrix::from_dense(&shifted, 1e-30),
                    krylov_dim,
                }
            }
        }
    }

    /// Apply the operator: compute A·v (Dense) or exp((G-r·I)·dt)·v (Krylov).
    ///
    /// Called at every timestep: J_k = operator.apply(J_{k+1} + dt·H3_{k+1})
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        match self {
            MatExpOperator::Dense { a_matrix } => dense_matvec(a_matrix, v),
            MatExpOperator::Sparse { g_shifted_dt, krylov_dim } => {
                let n = g_shifted_dt.nrows;
                let matvec = |x: &[f64]| g_shifted_dt.matvec(x);
                expm_krylov_action(&matvec, v, n, *krylov_dim, 1e-12)
            }
        }
    }

    /// Choose strategy based on problem size.
    pub fn suggest_strategy(nm: usize) -> MatExpStrategy {
        if nm <= 1500 {
            MatExpStrategy::DensePade
        } else {
            MatExpStrategy::Krylov { krylov_dim: 40.min(nm) }
        }
    }
}

fn dense_matvec(a: &DMatrix<f64>, x: &[f64]) -> Vec<f64> {
    debug_assert_eq!(x.len(), a.ncols());
    (a * DVector::from_column_slice(x)).as_slice().to_vec()
}

#[inline]
fn vec_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

#[inline]
fn vec_norm(v: &[f64]) -> f64 {
    vec_dot(v, v).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expm_zero_is_identity() {
        let n = 5;
        let result = expm(&DMatrix::zeros(n, n));
        let id = DMatrix::<f64>::identity(n, n);
        for i in 0..n {
            for j in 0..n {
                assert!((result[(i, j)] - id[(i, j)]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn expm_diagonal() {
        let diag = [1.0, -2.0, 0.5, 3.0_f64];
        let n = diag.len();
        let mut a = DMatrix::zeros(n, n);
        for i in 0..n {
            a[(i, i)] = diag[i];
        }
        let result = expm(&a);
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { diag[i].exp() } else { 0.0 };
                assert!((result[(i, j)] - expected).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn expm_nilpotent_2x2() {
        let mut a = DMatrix::zeros(2, 2);
        a[(0, 1)] = 1.0;
        let r = expm(&a);
        assert!((r[(0, 0)] - 1.0).abs() < 1e-14);
        assert!((r[(0, 1)] - 1.0).abs() < 1e-14);
        assert!((r[(1, 0)] - 0.0).abs() < 1e-14);
        assert!((r[(1, 1)] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn expm_generator_is_stochastic() {
        let mut g = DMatrix::zeros(3, 3);
        g[(0, 0)] = -2.0; g[(0, 1)] = 1.5; g[(0, 2)] = 0.5;
        g[(1, 0)] = 1.0;  g[(1, 1)] = -3.0; g[(1, 2)] = 2.0;
        g[(2, 0)] = 0.3;  g[(2, 1)] = 0.7;  g[(2, 2)] = -1.0;

        let p = expm(&(&g * 0.5));

        for i in 0..3 {
            let row_sum: f64 = (0..3).map(|j| p[(i, j)]).sum();
            assert!((row_sum - 1.0).abs() < 1e-12, "row {} sum = {}", i, row_sum);
            for j in 0..3 {
                assert!(p[(i, j)] >= -1e-14, "p[{},{}] = {}", i, j, p[(i, j)]);
            }
        }
    }

    #[test]
    fn krylov_matches_dense() {
        let n = 10;
        let mut g = DMatrix::zeros(n, n);
        for i in 0..n {
            let up = if i < n - 1 { 1.5 + 0.1 * i as f64 } else { 0.0 };
            let dn = if i > 0 { 1.0 + 0.05 * i as f64 } else { 0.0 };
            if i < n - 1 { g[(i, i + 1)] = up; }
            if i > 0 { g[(i, i - 1)] = dn; }
            g[(i, i)] = -(up + dn);
        }

        let a = &g * 0.3;
        let dense_a = expm(&a);
        let v: Vec<f64> = (0..n).map(|i| ((i as f64 + 1.0) as f64).sin()).collect();

        let dense_result = dense_matvec(&dense_a, &v);
        let krylov_result = expm_action_dense(&a, &v, n);

        let max_diff = dense_result
            .iter()
            .zip(&krylov_result)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        assert!(max_diff < 1e-10, "Krylov vs Dense max diff: {:.2e}", max_diff);
    }

    #[test]
    fn expm_large_norm_generator() {
        let mut g = DMatrix::zeros(3, 3);
        g[(0, 0)] = -100.0; g[(0, 1)] = 60.0;  g[(0, 2)] = 40.0;
        g[(1, 0)] = 30.0;   g[(1, 1)] = -80.0; g[(1, 2)] = 50.0;
        g[(2, 0)] = 20.0;   g[(2, 1)] = 40.0;  g[(2, 2)] = -60.0;

        let p = expm(&g);
        for i in 0..3 {
            let row_sum: f64 = (0..3).map(|j| p[(i, j)]).sum();
            assert!((row_sum - 1.0).abs() < 1e-10, "row {} sum = {}", i, row_sum);
        }
    }

    #[test]
    #[ignore]
    fn expm_ctmc_size_benchmark() {
        use std::time::Instant;
        for nm in [100, 500, 1000, 2500] {
            let mut g = DMatrix::zeros(nm, nm);
            for i in 0..nm {
                let up = if i < nm - 1 { 1.5 } else { 0.0 };
                let dn = if i > 0 { 1.0 } else { 0.0 };
                if i < nm - 1 { g[(i, i + 1)] = up; }
                if i > 0 { g[(i, i - 1)] = dn; }
                g[(i, i)] = -(up + dn);
            }
            let r = 0.05;
            let dt = 0.01;

            let t0 = Instant::now();
            let a = expm_generator(&g, r, dt);
            let t_dense = t0.elapsed();

            let v: Vec<f64> = (0..nm).map(|i| (i as f64 / nm as f64).max(0.0)).collect();
            let a_krylov = (&g - DMatrix::identity(nm, nm) * r) * dt;
            let t0 = Instant::now();
            let _ = expm_action_dense(&a_krylov, &v, 40.min(nm));
            let t_krylov = t0.elapsed();

            // Rows of exp((G-rI)*dt) should each sum to exp(-r*dt)
            let expected_sum = (-r * dt).exp();
            let max_err = (0..nm)
                .map(|i| ((0..nm).map(|j| a[(i, j)]).sum::<f64>() - expected_sum).abs())
                .fold(0.0_f64, f64::max);

            println!(
                "NM={nm:5}: Dense {t_dense:.2?}  Krylov {t_krylov:.2?}  row-sum err {max_err:.2e}"
            );
        }
    }
}
