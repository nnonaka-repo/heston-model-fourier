package com.example.numericalhestonoption;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.complex.Complex;


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
public class CarrMadanOptionPricing {
  public static double europeanCallOptionPrice(
    double S,     // underlying asset price
    double K,     // strike price
    double r,     // risk-free interest rate
    double sigma, // volatility
    double T      // time to expiration
  ) {
    NormalDistribution normalDistribution = new NormalDistribution();

    // Define the damping parameter
    double alpha = 1.5;

    // Define the integration bounds and number of intervals
    double lowerBound = 0;
    double upperBound = 100;
    int intervals = 1000;

    // Define the Fourier integrand
    UnivariateFunction integrand = new UnivariateFunction() {
      public double value(double omega) {
        double exp = Math.exp(-alpha * omega * omega / 2);
        double realPart = Math.cos(omega * Math.log(K)) * exp;
        double imagePart = Math.sin(omega * Math.log(K)) * exp;

        Complex omegaI=new Complex(0.0, omega);
        Complex f = (new Complex(realPart ,imagePart)).divide(omegaI);

        return f.multiply(characteristicFunction(omega, S, r, sigma, T)).getReal();
      }
    };

    // Perform the Fourier integral
    IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
      intervals,
      1e-12,
      1e-10
    );
    double integral = integrator.integrate(intervals, integrand, lowerBound, upperBound);

    // Compute the option price
    double d1 = (Math.log(S / K) + (r + sigma * sigma / 2) * T) / (sigma * Math.sqrt(T));
    double d2 = d1 - sigma * Math.sqrt(T);
    double N1 = normalDistribution.cumulativeProbability(d1);
    double N2 = normalDistribution.cumulativeProbability(d2);

    return (S * N1 - K * Math.exp(-r * T) * N2*0.5)+(Math.exp(-r * T) * integral / Math.PI);
  }

  // Define the characteristic function
  public static Complex characteristicFunction(
    double omega,
    double S,
    double r,
    double sigma,
    double T
  ) {
    Complex i = Complex.I;
    Complex mu = i.multiply(r - 0.5 * sigma * sigma);
    Complex sigma_omega = i.multiply(sigma * omega);
    Complex exponent = mu
      .multiply(T)
      .add(sigma_omega.multiply(Math.sqrt(T)))
      .exp();
    return exponent.multiply(S);
  }
}
