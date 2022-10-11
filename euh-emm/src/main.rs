// use plotly::Histogram;
// use plotly::Plot;
use rand::distributions::Distribution;
use statrs::distribution::{Binomial, Continuous, Normal};

fn one_mix(x_i: f64, mu_k: f64, sigma_k: f64, pi_k: f64) -> f64 {
    let norm = Normal::new(mu_k, sigma_k).unwrap();
    pi_k * norm.pdf(x_i)
}

fn log_lik<const K: usize>(x: &[f64], mu: &[f64; K], sigma: &[f64; K], pi: &[f64; K]) -> f64 {
    x.iter()
        .map(|&x_i| {
            mu.iter()
                .zip(sigma.iter())
                .zip(pi.iter())
                .map(|((&mu_k, &sigma_k), &pi_k)| {
                    pi_k * Normal::new(mu_k, sigma_k).unwrap().pdf(x_i)
                })
                .sum()
        })
        .map(f64::ln)
        .sum()
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
            x.iter()
                .map(|&x_i| x_i * gamma(x_i, k, mu, sigma, pi))
                .sum::<f64>()
                / gamma_sum[k]
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    let new_sigma: [f64; K] = (0..K)
        .map(|k| {
            x.iter()
                .map(|&x_i| (x_i - mu[k]).powi(2) * gamma(x_i, k, mu, sigma, pi))
                .sum::<f64>()
                / gamma_sum[k]
        })
        .map(f64::sqrt)
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
        [mu, sigma, pi] = m_step(x, &mu, &sigma, &pi);
    }
    [mu, sigma, pi]
}

fn main() {
    let num_samples = 5000;
    let num_steps = 100;
    let mut rng = rand::thread_rng();
    let binom = Binomial::new(0.75, 1).unwrap();
    let n1 = Normal::new(5.0, f64::sqrt(1.5)).unwrap();
    let n2 = Normal::new(10.0, f64::sqrt(2.0)).unwrap();

    // Generate samples
    let x: Vec<_> = (0..num_samples)
        .map(|_| {
            let z = binom.sample(&mut rng);
            if z <= 0.5 {
                n1.sample(&mut rng)
            } else {
                n2.sample(&mut rng)
            }
        })
        .collect();

    // let trace = Histogram::new(x.clone()).name("h");
    // let mut plot = Plot::new();
    // plot.add_trace(trace);
    // plot.show();

    // println!("{x:#?}");

    // Initial parameters
    let mu_0 = [-1.0, 20.0];
    let sigma_0 = [f64::sqrt(1.7), f64::sqrt(1.7)];
    let pi_0 = [0.5, 0.5];

    let init_log_lik = log_lik(&x, &mu_0, &sigma_0, &pi_0);
    // Run the algo
    let begin = std::time::Instant::now();
    let [mu, sigma, pi] = em_algo(&x, &mu_0, &sigma_0, &pi_0, num_steps);
    let elapsed = begin.elapsed();

    let final_log_lik = log_lik(&x, &mu, &sigma, &pi);

    println!("Took {} seconds", elapsed.as_secs_f32());
    println!("Initial log likelihood: {init_log_lik:.3}");
    println!("After {num_steps} steps, final log likelihood: {final_log_lik:.3}");
    println!("mu: {mu:.3?}, sigma: {sigma:.3?}, pi: {pi:.3?}");

    // Gibbs sampling on multivariate normal
    // let mvn = MultivariateNormal::new(vec![1.0, 1.0], vec![1.0, 4.0, 1.0, 4.0]).unwrap();
}
