//! This file contains implementations for the learning of a graphical bayesian graph using
//! conditional expectation maximisation
#![warn(
    rust_2018_idioms,
    rust_2021_compatibility,
    missing_debug_implementations,
    missing_docs
)]
#![feature(generic_const_exprs)]

// use itertools::Itertools;
use na::{Cholesky, Const, DimMin, SMatrix, SVector, Scalar};
use nalgebra as na;
use rand::{distributions::Distribution, Rng};
use statrs::distribution::{Continuous, Normal};
use thiserror::Error;

type Matrix<const D: usize> = SMatrix<f64, D, D>;

struct EStep<const D: usize> {
    v0: f64,
    v1: f64,
    n0: Normal,
    n1: Normal,
}

impl<const D: usize> EStep<D> {
    fn new(v0: f64, v1: f64) -> Self {
        Self {
            v0,
            v1,
            n0: Normal::new(0.0, v0).unwrap(),
            n1: Normal::new(0.0, v1).unwrap(),
        }
    }

    fn compute(&self, a: &Matrix<D>, pi: f64) -> (Matrix<D>, Matrix<D>) {
        let p: Matrix<D> = SMatrix::from_iterator(a.iter().map(|&a_jk| {
            let alpha = self.n1.pdf(a_jk) * pi;
            let beta = self.n0.pdf(a_jk) * (1.0 - pi);
            alpha / (alpha + beta)
        }));
        let d: Matrix<D> = p.map(|p_jk| (1.0 - p_jk) / self.v0.powi(2) + p_jk / self.v1.powi(2));

        (p, d)
    }
}

struct CMStep<const D: usize> {
    alpha: f64,
    beta: f64,
    lambda: f64,
    s: Matrix<D>,
    v0: f64,
    v1: f64,
    n: usize,
}

impl<const D: usize> CMStep<D> {
    fn new<const N: usize>(
        alpha: f64,
        beta: f64,
        lambda: f64,
        v0: f64,
        v1: f64,
        x: &SMatrix<f64, N, D>,
    ) -> Self {
        let s = x.tr_mul(x);
        CMStep {
            s,
            alpha,
            beta,
            lambda,
            v0,
            v1,
            n: N,
        }
    }

    fn compute(&self, a: &Matrix<D>, p: &Matrix<D>, d: &mut Matrix<D>) -> (Matrix<D>, f64)
    where
        [(); D - 1]:,
    {
        let mut s = self.s.clone_owned();
        let sum_delta = p.upper_triangle().sum();
        let pi = (self.alpha + sum_delta - 1.0)
            / (self.alpha + self.beta + (D * (D - 1) / 2) as f64 - 2.0);
        let mut omega = a.clone_owned();

        for t in 0..D {
            // Swap columns and calculate components based off original omega
            omega.swap_columns(t, D - 1);
            omega.swap_rows(t, D - 1);
            d.swap_columns(t, D - 1);
            d.swap_rows(t, D - 1);
            s.swap_columns(t, D - 1);
            s.swap_rows(t, D - 1);

            // Get relevant values
            let s_22 = s.get((D - 1, D - 1)).unwrap();
            let s_12 = s.slice((0, D - 1), (D - 1, 1));
            let d_12 = d.fixed_slice(0, D - 1);

            let omega_11_inv = omega.slice((0, 0), (D - 1, D - 1)).try_inverse().unwrap();
            let c: Matrix<{ D - 1 }> =
                (s_22 + self.lambda) * omega_11_inv.clone() + Matrix::from_diagonal(&d_12);

            // Update column upper part
            let omega_12_lp1 = -c.try_inverse().unwrap() * s_12;
            omega
                .slice_mut((0, D - 1), (D - 1, 1))
                .copy_from(&omega_12_lp1);
            omega
                .slice_mut((D - 1, 0), (1, D - 1))
                .copy_from(&omega_12_lp1.transpose());

            // Update column diagonal bit
            let omega_22_lp1 = (omega_12_lp1.tr_mul(&omega_11_inv) * omega_12_lp1)
                .add_scalar(self.n as f64 / (self.lambda + s_22));
            omega
                .slice_mut((D - 1, D - 1), (1, 1))
                .copy_from(&omega_22_lp1);

            // Swap back column, so that we are ready for next iteration
            omega.swap_columns(t, D - 1);
            omega.swap_rows(t, D - 1);
            d.swap_columns(t, D - 1);
            d.swap_rows(t, D - 1);
            s.swap_columns(t, D - 1);
            s.swap_rows(t, D - 1);
        }
        (omega, pi)
    }

    fn log_likelihood(&self, omega: &Matrix<D>, pi: f64, p: &Matrix<D>, d: &Matrix<D>) -> f64
    where
        Const<D>: DimMin<Const<D>, Output = Const<D>>,
    {
        let det = omega.determinant();
        let sns_bit = {
            let mut mat = omega.zip_map(d, |o_ij, d_ij| d_ij * o_ij.powi(2));
            mat.fill_lower_triangle(0.0, 0);
            -0.5 * (mat.sum() + self.lambda * omega.diagonal().sum())
        };
        let ber_bit = {
            let mut mat = p.map(|p_ij| (1.0 - pi).ln() + p_ij * (pi / (1.0 - pi)).ln());
            mat.fill_lower_triangle(0.0, 0);
            mat.sum()
        };
        let misc_bit = (self.alpha - 1.0) * pi.ln()
            + (self.beta - 1.0) * dbg!((1.0 - pi).ln())
            + 0.5 * self.n as f64 * dbg!(det.ln())
            - 0.5 * dbg!((self.s * omega).trace());

        dbg!(misc_bit) + dbg!(ber_bit) - dbg!(sns_bit)
    }
}

struct MultivariateNormal<const D: usize> {
    mean: SVector<f64, D>,
    chol_cov: Cholesky<f64, Const<D>>,
}

impl<const D: usize> MultivariateNormal<D> {
    fn new(mean: Vec<f64>, prec: Vec<f64>) -> Self {
        let mean = SVector::from_vec(mean);
        let prec = Matrix::<D>::from_vec(prec);
        let cov = prec.try_inverse().unwrap();
        let chol_cov = cov.cholesky().unwrap();
        MultivariateNormal { mean, chol_cov }
    }
}

impl<const D: usize> Distribution<SVector<f64, D>> for MultivariateNormal<D> {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> SVector<f64, D> {
        let norm = Normal::new(0.0, 1.0).unwrap();
        let z = SVector::from_iterator((0..D).map(|_| norm.sample(rng)));
        self.chol_cov.l() * z + self.mean
    }
}

const M_SIZE: usize = 20;
const SAMPLE_SIZE: usize = 100;
const NUM_STEPS: usize = 50;

fn main() {
    let mut ar2 = Matrix::<M_SIZE>::zeros();
    ar2.fill_with_identity();
    ar2.fill_upper_triangle(0.5, 1);
    ar2.fill_upper_triangle(0.25, 2);
    ar2.fill_upper_triangle(0.0, 3);
    ar2.fill_lower_triangle(0.5, 1);
    ar2.fill_lower_triangle(0.25, 2);
    ar2.fill_lower_triangle(0.0, 3);

    let prec: Vec<f64> = ar2.iter().cloned().collect();
    let mut rng = rand::thread_rng();
    let mvn = MultivariateNormal::<M_SIZE>::new(vec![0.0; M_SIZE], prec);
    let x: SMatrix<f64, SAMPLE_SIZE, M_SIZE> = mvn
        .sample_iter(&mut rng)
        .take(SAMPLE_SIZE)
        .enumerate()
        .fold(
            SMatrix::<f64, SAMPLE_SIZE, M_SIZE>::zeros(),
            |mut mat, (i, x_i)| {
                // println!("{:?}", x_i.shape());
                // println!("{x_i}");
                mat.set_row(i, &x_i.transpose());
                mat
            },
        );

    let alpha = 1.0;
    let beta = 1.0;
    let lambda = 1.0;
    let mut pi = 0.5;
    let v0 = 0.1;
    let v1 = 10.0;
    let mut omega = Matrix::from_diagonal_element(4.0);
    let e_step = EStep::new(v0, v1);
    let cm_step = CMStep::new(alpha, beta, lambda, v0, v1, &x);

    for _ in 0..NUM_STEPS {
        let (p, mut d) = e_step.compute(&omega, pi);
        (omega, pi) = cm_step.compute(&omega, &p, &mut d);
        println!("{}", cm_step.log_likelihood(&omega, pi, &p, &d));
    }
    println!("{pi}");
    // println!("{omega:#?}");
    println!("{:?}", e_step.compute(&omega, pi).0);
}

/// An error enum for the methods in this file
#[derive(Debug, Error)]
pub enum SNSError {
    /// This error comes directly from `statrs::StatsError`
    #[error(transparent)]
    StatsError(#[from] statrs::StatsError),
}
