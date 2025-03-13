package com.nnonaka.option.price.service;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;

import org.apache.commons.math3.exception.MaxCountExceededException;
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
public class CarrMadanOptionPricingService {
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

    // Define the integration bounds and number of intervals
    double lowerBound = 0;
    double upperBound = 100;
    int intervals = 1000; // Small number of intervals
    int maxEvaluations = intervals * 10; // Adjust based on the number of intervals


    // Define the Fourier integrand
    UnivariateFunction integrand = new UnivariateFunction() {
      public double value(double omega) {
        Complex complexOmega = new Complex(0, omega);
        Complex characteristic = characteristicFunction(complexOmega, S, r, sigma, T);
        Complex dampedChar = complexOmega.multiply(-alpha * omega / 2.0).exp().multiply(characteristic);
        Complex f = dampedChar.divide(complexOmega);
        return f.getReal();
      }
    };

    // Perform the Fourier integral
    IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
            intervals,
            1e-10,
            1e-10);
    double integral;
    try {
      integral = integrator.integrate(maxEvaluations, integrand, lowerBound, upperBound);
    } catch (MaxCountExceededException e) {
      // Handle the exception here
      // You can choose to return a default value or take any other appropriate action
      System.out.printf("Integration error: %s.", e.getMessage());
      return String.format("Integration error: %s.", e.getMessage());
    }



    // Compute the option price
    double d1 = (Math.log(S / K) + (r + sigma * sigma / 2) * T) / (sigma * Math.sqrt(T));
    double d2 = d1 - sigma * Math.sqrt(T);
    double N1 = normalDistribution.cumulativeProbability(d1);
    double N2 = normalDistribution.cumulativeProbability(d2);
    double callPrice = (S * N1 - K * Math.exp(-r * T) * N2*0.5) + (Math.exp(-r * T) * integral / Math.PI);
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
