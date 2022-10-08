use rand::distributions::Distribution;
use statrs::distribution::{Binomial, Continuous, Normal};

fn one_mix(x_i: f64, mu_k: f64, sigma_k: f64, pi_k: f64) -> f64 {
    let norm = Normal::new(mu_k, sigma_k.powi(2)).unwrap();
    pi_k * norm.pdf(x_i)
}

fn gamma(x_i: f64, k: usize, mu: &[f64], sigma: &[f64], pi: &[f64]) -> f64 {
    let denominator: f64 = mu
        .iter()
        .zip(sigma.iter())
        .zip(pi.iter())
        .map(|((mu_k, sigma_k), pi_k)| one_mix(x_i, *mu_k, *sigma_k, *pi_k))
        .sum();
    one_mix(x_i, mu[k], sigma[k], pi[k]) / denominator
}

fn m_step<const K: usize>(
    x: &[f64],
    mu: &[f64; K],
    sigma: &[f64; K],
    pi: &[f64; K],
) -> [[f64; K]; 3] {
    let gamma_sum: Vec<f64> = (0..K)
        .map(|k| x.iter().map(|x_i| gamma(*x_i, k, mu, sigma, pi)).sum())
        .collect();

    let new_mu: [f64; K] = (0..K)
        .map(|k| {
            1.0 / gamma_sum[k]
                * x.iter()
                    .map(|&x_i| x_i * gamma(x_i, k, mu, sigma, pi))
                    .sum::<f64>()
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    let new_sigma: [f64; K] = (0..K)
        .map(|k| {
            1.0 / gamma_sum[k]
                * x.iter()
                    .map(|&x_i| (x_i - mu[k]).powi(2) * gamma(x_i, k, mu, sigma, pi))
                    .sum::<f64>()
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let new_pi: [f64; K] = (0..K)
        .map(|k| gamma_sum[k] / x.len() as f64)
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    [new_mu, new_sigma, new_pi]
}

fn em_algo<const K: usize>(
    x: &[f64],
    mu: &[f64; K],
    sigma: &[f64; K],
    pi: &[f64; K],
    num_steps: u32,
) -> [[f64; K]; 3] {
    let mut mu = *mu;
    let mut sigma = *sigma;
    let mut pi = *pi;
    for _ in 0..num_steps {
        [mu, sigma, pi] = dbg!(m_step(x, &mu, &sigma, &pi));
    }
    [mu, sigma, pi]
}

fn main() {
    let num_samples = 5000;
    let mut rng = rand::thread_rng();
    let binom = Binomial::new(0.75, 1).unwrap();
    let x: Vec<_> = (0..num_samples)
        .map(|_| {
            let z = binom.sample(&mut rng);
            let mean = if z <= 0.5 { 40.0 } else { 64.0 };
            Normal::new(mean, 2.5).unwrap().sample(&mut rng)
        })
        .collect();
    // println!("prices: {prices:#?}");

    let num_steps = 4;
    let mu_0 = [60.0, 67.0];
    let sigma_0 = [1.0, 1.1];
    let pi_0 = [0.6, 0.4];

    println!(
        "after {num_steps} ... {:?}",
        em_algo(&x, &mu_0, &sigma_0, &pi_0, num_steps)
    );
}
