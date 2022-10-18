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
use na::{Cholesky, Const, SMatrix, SVector};
use nalgebra as na;
use rand::{distributions::Distribution, Rng};
use statrs::distribution::{Continuous, Normal};
use thiserror::Error;

type Matrix<const D: usize> = SMatrix<f64, D, D>;
type BMatrix<const D: usize> = SMatrix<bool, D, D>;

// struct MatrixPrior {
//     n0: Normal,
//     n1: Normal,
//     exp: Exp,
//     pi: f64,
// }
//
// impl MatrixPrior {
//     fn new(v0: f64, v1: f64, pi: f64, lambda: f64) -> Result<Self, SNSError> {
//         Ok(Self {
//             n0: Normal::new(0.0, v0)?,
//             n1: Normal::new(0.0, v1)?,
//             exp: Exp::new(lambda / 2.0)?,
//             pi,
//         })
//     }
// }
//
// impl<const D: usize> Continuous<Matrix<D>, f64> for MatrixPrior {
//     fn pdf(&self, a: Matrix<D>) -> f64 {
//         if a.cholesky().is_none() {
//             // check if `a` is positive definite
//             return 0.0;
//         }
//
//         // Get elements of strictly upper triangular part of `a`
//         let sns_prod: f64 = (0..D)
//             .cartesian_product(0..D)
//             .filter(|(i, j)| i < j)
//             .filter_map(|(i, j)| a.get((i, j)))
//             .map(|&a_ij| (1.0 - self.pi) * self.n0.pdf(a_ij) + self.pi * self.n1.pdf(a_ij))
//             .product();
//         // Compute exponentially distributed part
//         let exp_prod = a.diagonal().map(|a_ii| self.exp.pdf(a_ii)).product();
//
//         sns_prod * exp_prod
//     }
//     fn ln_pdf(&self, x: Matrix<D>) -> f64 {
//         // TODO: Rewrite it as the whole logged formula
//         self.pdf(x).ln()
//     }
// }
//
// struct MatrixGivenGraph<const D: usize> {
//     z: BMatrix<D>,
//     n0: Normal,
//     n1: Normal,
//     exp: Exp,
// }
//
// impl<const D: usize> Continuous<Matrix<D>, f64> for MatrixGivenGraph<D> {
//     fn pdf(&self, a: Matrix<D>) -> f64 {
//         let normal_prod: f64 = (0..D)
//             .cartesian_product(0..D)
//             .filter(|(i, j)| i < j)
//             .filter_map(|(i, j)| a.get((i, j)).zip(self.z.get((i, j))))
//             .map(|(&a_ij, &z_ij)| {
//                 if z_ij {
//                     self.n1.pdf(a_ij)
//                 } else {
//                     self.n0.pdf(a_ij)
//                 }
//             })
//             .product();
//         let exp_prod = a.diagonal().map(|a_ii| self.exp.pdf(a_ii)).product();
//
//         normal_prod * exp_prod
//     }
//
//     fn ln_pdf(&self, a: Matrix<D>) -> f64 {
//         self.pdf(a).ln()
//     }
// }
//
// struct GraphPrior {
//     pi: f64,
// }
//
// impl<const D: usize> Continuous<BMatrix<D>, f64> for GraphPrior {
//     fn pdf(&self, z: BMatrix<D>) -> f64 {
//         (0..D)
//             .cartesian_product(0..D)
//             .filter(|(i, j)| i < j)
//             .filter_map(|(i, j)| z.get((i, j)))
//             .map(|&z_ij| self.pi.powi(z_ij as i32) * (1.0 - self.pi).powi(z_ij as i32))
//             .product()
//     }
//     fn ln_pdf(&self, z: BMatrix<D>) -> f64 {
//         self.pdf(z).ln()
//     }
// }

struct EStep<const D: usize> {
    // p: Matrix<D>,
    d: Matrix<D>,
}

impl<const D: usize> EStep<D> {
    fn compute(a: &Matrix<D>, pi: f64, v0: f64, v1: f64) -> Self {
        let p: Matrix<D> = SMatrix::from_iterator(a.iter().map(|&a_jk| {
            let alpha = Normal::new(0.0, v1).unwrap().pdf(a_jk) * pi;
            let beta = Normal::new(0.0, v0).unwrap().pdf(a_jk) * (1.0 - pi);
            alpha / (alpha + beta)
        }));
        let d: Matrix<D> = p.map(|p_jk| (1.0 - p_jk) / v0.powi(2) + p_jk / v1.powi(2));

        Self { d }
    }
}

struct CMStep {
    alpha: f64,
    beta: f64,
    lambda: f64,
    pi: f64,
    v0: f64,
    v1: f64,
}

impl CMStep {
    fn new(alpha: f64, beta: f64, lambda: f64, pi: f64, v0: f64, v1: f64) -> Self {
        CMStep {
            alpha,
            beta,
            lambda,
            pi,
            v0,
            v1,
        }
    }

    fn compute<const D: usize, const N: usize>(
        &self,
        x: &SMatrix<f64, N, D>,
        z: &BMatrix<D>,
        a: &Matrix<D>,
    ) -> (f64, Matrix<D>)
    where
        [(); D - 1]:,
    {
        let p = D;
        let n = N as f64;
        let s_z = z.map(|z_ij| z_ij as i32).sum() as f64;
        let new_pi =
            (self.alpha + s_z - 1.0) / (self.alpha + self.beta + (p * (p - 1) / 2) as f64 - 2.0);
        let s = x.tr_mul(x);
        let d = EStep::compute(a, self.pi, self.v0, self.v1).d;

        println!("d: {d}");

        let mut omega = a.clone_owned();
        let s_22 = s.get((D - 1, D - 1)).unwrap();
        let s_12 = s.index((0..D - 1, D - 1));
        let d_12_inv = Matrix::<{ D - 1 }>::from_diagonal(&d.fixed_slice(0, D - 1))
            .try_inverse()
            .unwrap();
        let mut new_omega = Matrix::<D>::zeros();

        for t in 0..D {
            omega.swap_columns(t, D - 1);
            let omega_11_inv = omega.slice((0, 0), (D - 1, D - 1)).try_inverse().unwrap();
            let c: Matrix<{ D - 1 }> = (s_22 + self.lambda) * omega_11_inv.clone() + d_12_inv;

            // update column upper part
            let omega_12_lp1 = c.try_inverse().unwrap() * s_12;
            new_omega
                .slice_mut((0, t), (D - 1, 1))
                .copy_from(&omega_12_lp1);

            // update column lower part
            let omega_22_lp1 = (omega_12_lp1.tr_mul(&omega_11_inv) * omega_12_lp1)
                .add_scalar(n / (self.lambda + s_22));
            new_omega
                .slice_mut((D - 1, t), (1, 1))
                .copy_from(&omega_22_lp1);
            omega.swap_columns(t, D - 1);
        }
        (new_pi, new_omega)
    }
}

struct MultivariateNormal<const D: usize> {
    // cov: Matrix<D>,
    // prec: Matrix<D>,
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

const M_SIZE: usize = 5;
const SAMPLE_SIZE: usize = 100;
const NUM_STEPS: usize = 4;

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

    let z = ar2.map(|a_ij| rng.gen_bool(a_ij));
    // dbg!(&x);
    // dbg!(&z);

    let alpha = 1.0;
    let beta = 1.0;
    let lambda = 1.0;
    let mut pi = 0.5;
    let v0 = 1.0;
    let v1 = 100.0;

    for _ in 0..NUM_STEPS {
        // println!("{omega}");
        (pi, ar2) = CMStep::new(alpha, beta, lambda, pi, v0, v1).compute(&x, &z, &ar2);
        println!("pi: {pi}, omega: {ar2}");
    }
    println!("{ar2}");
}

/// An error enum for the methods in this file
#[derive(Debug, Error)]
pub enum SNSError {
    /// This error comes directly from `statrs::StatsError`
    #[error(transparent)]
    StatsError(#[from] statrs::StatsError),
}
