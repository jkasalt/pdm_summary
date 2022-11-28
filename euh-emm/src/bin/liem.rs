//! This file contains implementations for the learning of a graphical bayesian graph using
//! conditional expectation maximisation
#![feature(generic_const_exprs)]

use approx::relative_eq;
use euh_emm::distributions::MultivariateNormal;
use na::{Const, DMatrix, DimMin, SMatrix};
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
        let p: Matrix<D> = SMatrix::from_row_iterator(a.iter().map(|&a_jk| {
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
    n: usize,
}

impl<const D: usize> CMStep<D> {
    fn new<const N: usize>(alpha: f64, beta: f64, lambda: f64, x: &SMatrix<f64, N, D>) -> Self {
        let s = x.tr_mul(x);
        CMStep {
            s,
            alpha,
            beta,
            lambda,
            n: N,
        }
    }

    fn compute(&self, a: &Matrix<D>, p: &Matrix<D>, d: &mut Matrix<D>) -> (Matrix<D>, f64)
    where
        [(); D - 1]:,
    {
        let mut s = self.s.clone_owned();
        let sum_delta = p.upper_triangle().sum() - p.diagonal().sum();
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
            let omega_22_lp1 = (&omega_12_lp1.tr_mul(&omega_11_inv) * omega_12_lp1)
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

        assert!(
            omega.cholesky().is_some(),
            "CMStep has failed to produce a positive definite matrix"
        );
        (omega, pi)
    }

    fn log_likelihood(&self, omega: &Matrix<D>, pi: f64, p: &Matrix<D>, d: &Matrix<D>) -> f64
    where
        Const<D>: DimMin<Const<D>, Output = Const<D>>,
    {
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
        let det = omega.determinant();
        let misc_bit = (self.alpha - 1.0) * pi.ln()
            + (self.beta - 1.0) * (1.0 - pi).ln()
            + 0.5 * self.n as f64 * det.ln()
            - 0.5 * (self.s * omega).trace();

        misc_bit + ber_bit + sns_bit
    }
}

/// Computes the true positive rate for an estimated posterior inclusion matrix
fn tpr<const D: usize>(p: &Matrix<D>, truth: &Matrix<D>, cutoff: f64) -> f64 {
    let total = p.map(|pij| (pij > cutoff) as i32).sum() as f64;
    let true_pos = p
        .zip_map(truth, |pij, tij| ((pij > cutoff) && (tij > 0.0)) as i32)
        .sum() as f64;

    dbg!(true_pos) / dbg!(total)
}

fn fpr<const D: usize>(p: &Matrix<D>, truth: &Matrix<D>, cutoff: f64) -> f64 {
    let total = p.map(|pij| (pij > cutoff) as i32).sum() as f64;
    let false_pos = p
        .zip_map(truth, |pij, tij| {
            ((pij > cutoff) && (relative_eq!(tij, 0.0))) as i32
        })
        .sum() as f64;

    false_pos / total
}

/// `rate` is the percentage of non-zero entries, up to some random variance
// fn prior_precision<const D: usize>(rng: &mut impl Rng) -> Matrix<D> {
//     assert!(
//         (0.0..=1.0).contains(&rate),
//         "rate has to be in the range of [0,1]"
//     );
//     let mat = Matrix::<D>::new_random();
//     let pos_def =
//         (dbg!(mat) * mat.transpose()).map(|aij| if rng.gen_bool(rate) { aij } else { 0.0 });
//     assert!(
//         pos_def.cholesky().is_some(),
//         "random_sparse_pos_def has failed to generate a positive definite matrix"
//     );
//     pos_def
// }

const M_SIZE: usize = 500;
const SAMPLE_SIZE: usize = 2000;
const NUM_STEPS: usize = 30;

fn main() {
    let file = std::fs::File::open("../r-bmg/x.csv").unwrap();
    let mut rdr = csv::Reader::from_reader(file);
    let x: Vec<f64> = rdr.records().fold(Vec::new(), |mut vec, line| {
        let line = line.unwrap();
        vec.extend(line.into_iter().filter_map(|l| l.parse::<f64>().ok()));
        vec
    });
    println!("{}", x.len());
    let x = SMatrix::<f64, SAMPLE_SIZE, M_SIZE>::from_row_slice(&x);
    let alpha = 1.0;
    let beta = 1.0;
    let lambda = 1.0;
    let mut pi = 0.5;
    let v0 = 0.1;
    let v1 = 10.0;
    let mut omega = Matrix::identity();
    let e_step = EStep::new(v0, v1);
    let cm_step = CMStep::new(alpha, beta, lambda, &x);

    for _ in 0..NUM_STEPS {
        let (p, mut d) = e_step.compute(&omega, pi);
        (omega, pi) = cm_step.compute(&omega, &p, &mut d);
        // println!("{}", cm_step.log_likelihood(&omega, pi, &p, &d));
    }
    // let (p, d) = e_step.compute(&omega, pi);
    println!("pi: {pi}, omega first row: {}", omega.row(0));
    // println!("p first row: {}", p.row(0));
    // println!("{omega:#?}");
    // let (p, _d) = e_step.compute(&omega, pi);
    // println!("tpr: {}, fpr: {}", tpr(&p, &ar2, 0.05), fpr(&p, &ar2, 0.5));
}

/// An error enum for the methods in this file
#[derive(Debug, Error)]
pub enum SNSError {
    /// This error comes directly from `statrs::StatsError`
    #[error(transparent)]
    StatsError(#[from] statrs::StatsError),
}

#[cfg(test)]
mod tests {
    use super::*;
}
