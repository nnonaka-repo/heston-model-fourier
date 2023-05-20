package com.nnonaka.option.price.service;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.distribution.NormalDistribution;

import org.springframework.stereotype.Service;

@Service
public class GarmanKohlhagenOptionFFTPricing {
    public String europeanCallOptionPrice(double spotPrice, double strikePrice, double riskFreeRate, double volatility, double timeToExpiration) {

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
            Complex characteristic = characteristicFunction(omega, spotPrice, riskFreeRate, volatility, timeToExpiration);
            phi[i] = characteristic;
        }

        // Perform the FFT on the characteristic function values
        FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] fft = transformer.transform(phi, TransformType.FORWARD);

        // Compute the option price using the FFT result
        double callPrice = 0.0;
        double aux = 0.0;
        for (int i = 0; i < N; i++) {
            Complex omega = new Complex(0, x[i]);
            Complex dampedChar = omega.multiply(-riskFreeRate * x[i] / 2.0).exp().multiply(fft[i]);
            Complex f = dampedChar.divide(omega);
            aux=f.getReal() * delta;
            if(!Double.isNaN(aux))
                callPrice += aux;
        }
        callPrice *= Math.exp(-riskFreeRate * timeToExpiration) / Math.PI;

        // Compute the remaining components of the option price formula
        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + volatility * volatility / 2) * timeToExpiration) / (volatility * Math.sqrt(timeToExpiration));
        double d2 = d1 - volatility * Math.sqrt(timeToExpiration);
        double N1 = normalCDF(d1);
        double N2 = normalCDF(d2);

        callPrice += spotPrice * N1 - strikePrice * Math.exp(-riskFreeRate * timeToExpiration) * N2;
        return  String.format("callPrice : %f", callPrice);
    }

    /*public Complex characteristicFunction(Complex omega, double spotPrice, double riskFreeRate, double volatility, double timeToExpiration) {
        Complex iOmega = omega.multiply(-1);
        Complex d1 = iOmega.add(0.5 * volatility * volatility * timeToExpiration).sqrt().add(iOmega.multiply(volatility * Math.sqrt(timeToExpiration)));
        Complex d2 = iOmega.subtract(0.5 * volatility * volatility * timeToExpiration).sqrt().subtract(iOmega.multiply(volatility * Math.sqrt(timeToExpiration)));
        Complex numerator = d1.exp().subtract(1).multiply(iOmega.multiply(2)).add(1).multiply(spotPrice * Math.exp(-riskFreeRate * timeToExpiration));
        Complex denominator = iOmega.multiply(iOmega.add(1)).multiply(d2).multiply(d2.exp());
        return numerator.divide(denominator);
    }*/
    public static Complex characteristicFunction(Complex omega, double S, double r, double sigma, double T) {
        Complex i = Complex.I;
        Complex mu = i.multiply(r - 0.5 * sigma * sigma);
        Complex sigmaOmega = i.multiply(sigma).multiply(omega);
        Complex exponent = mu.multiply(T).add(sigmaOmega.multiply(Math.sqrt(T))).exp();
        return exponent.multiply(S);
    }
    public double normalCDF(double x) {
        return 0.5 * (1 + Erf.erf(x / Math.sqrt(2)));
    }
}
