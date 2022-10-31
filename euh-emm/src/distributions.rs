use nalgebra::{SVector, SMatrix, Cholesky, Const};
use statrs::distribution::Normal;
use rand::{Rng, distributions::Distribution};
type Matrix<const D: usize> = SMatrix<f64, D, D>;

pub struct MultivariateNormal<const D: usize> {
    mean: SVector<f64, D>,
    chol_cov: Cholesky<f64, Const<D>>,
}

impl<const D: usize> MultivariateNormal<D> {
    pub fn new(mean: Vec<f64>, prec: Vec<f64>) -> Self {
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

pub struct Wishart<const D: usize> {
    d: Matrix<D>,
    delta: f64,
}


