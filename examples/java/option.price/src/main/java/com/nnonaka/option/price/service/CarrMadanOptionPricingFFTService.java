package com.nnonaka.option.price.service;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import org.springframework.stereotype.Service;


//C(K,T) = e^(-rT) * (1/2 + (1/π) * ∫[0, ∞] Re[(e^(-iωlog(K)) / (iω)) * F(ω, alpha) ]dω)
/*
 * This code defines a function europeanCallOptionPrice
 * and returns the price of a European call option using the Carr-Madan formula.
 *
 * The function uses the characteristicFunction function to compute the damped characteristic function,
 * and performs a numerical Fourier integral using the Apache Commons Math library.
 *
 * The resulting option price is returned as a double value.
 */
@Service
public class CarrMadanOptionPricingFFTService {
  public String europeanCallOptionPrice(
    double K,     // strike price
    double sigma, // volatility
    double S,     // underlying asset price
    double r,     // risk-free interest rate
    double T      // time to expiration
  ) {
    NormalDistribution normalDistribution = new NormalDistribution();

    // Define the damping parameter
    double alpha = 1.5;

    int N = 4096; // Number of discretization points
    double delta = 0.01; // Discretization step size
    double[] x = new double[N]; // Discretization points

    // Compute the discretization points using the FFT frequencies
    for (int i = 0; i < N; i++) {
      x[i] = i * delta - (N / 2) * delta;
    }

    // Compute the characteristic function values at the discretization points
    Complex[] phi = new Complex[N];
    for (int i = 0; i < N; i++) {
      Complex omega = new Complex(0, x[i]);
      Complex characteristic = characteristicFunction(omega, S, r, sigma, T);
      phi[i] = characteristic;
    }

    // Perform the FFT on the characteristic function values
    FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
    Complex[] fft = transformer.transform(phi, TransformType.FORWARD);

    // Compute the integrand values
    double[] integrand = new double[N];
    for (int i = 0; i < N; i++) {
      Complex omega = new Complex(0, x[i]);
      Complex dampedChar = omega.multiply(-alpha * x[i] / 2.0).exp().multiply(fft[i]);
      Complex f = dampedChar.divide(omega);
      integrand[i] = f.getReal();
    }

    // Compute the integral using the trapezoidal rule
    double integral = 0.0;
    double aux = 0.0;
    for (int i = 0; i < N - 1; i++) {
      aux=(integrand[i] + integrand[i + 1]) * delta / 2.0;
      if(!Double.isNaN(aux))
        integral += aux;
    }

    // Compute the remaining components of the option price formula
    double d1 = (Math.log(S / K) + (r + sigma * sigma / 2) * T) / (sigma * Math.sqrt(T));
    double d2 = d1 - sigma * Math.sqrt(T);
    double N1 = normalDistribution.cumulativeProbability(d1);
    double N2 = normalDistribution.cumulativeProbability(d2);
    double callPrice = S * N1 - K * Math.exp(-r * T) * N2;
    callPrice += Math.exp(-r * T) * integral / Math.PI;

    return  String.format("callPrice : %f", callPrice);
  }

  // Define the characteristic function
  public static Complex characteristicFunction(Complex omega, double S, double r, double sigma, double T) {
    Complex i = Complex.I;
    Complex mu = i.multiply(r - 0.5 * sigma * sigma);
    Complex sigmaOmega = i.multiply(sigma).multiply(omega);
    Complex exponent = mu.multiply(T).add(sigmaOmega.multiply(Math.sqrt(T))).exp();
    return exponent.multiply(S);
  }


}
