package com.nnonaka.option.price.service;

import java.lang.Math;

import org.springframework.stereotype.Service;

// Updated Heston Parameters for February 21, 2025

// 1. S0 - Spot Price (USD/BRL)
// - Original: 5.80
// - Today: USD/BRL forecasts for 2025 range from 5.7 to 6.0 (e.g., 5.769 today per WalletInvestor, 6.5497 by year-end per Traders Union). Since it’s February, we’ll use a mid-range value.
// - New Value: S0 = 5.75 (slightly stronger BRL than your original).

// 2. K - Strike Price
// - Original: 5.85
// - This is fixed for the option you’re pricing.
// - New Value: K = 5.85 (no change).

// 3. T - Time to Maturity
// - Original: 1.0 (1 year)
// - Still 1 year from today (Feb 21, 2025 to Feb 21, 2026).
// - New Value: T = 1.0 (no change).

// 4. r - USD Risk-Free Rate
// - Original: 0.045 (4.5%)
// - Today: 1-year US Treasury yield in Feb 2025 is likely 4.0-4.5% (post-2024 Fed adjustments). We’ll use a slight dip from your original.
// - New Value: r = 0.0425 (4.25%).

// 5. q - BRL Risk-Free Rate
// - Original: 0.105 (10.5%)
// - Today: Brazil’s SELIC rate is around 10-11% in Feb 2025 (inflation moderating but high). Slightly up from your original.
// - New Value: q = 0.1075 (10.75%).

// 6. v0 - Initial Variance
// - Original: 0.04 (volatility = square root of 0.04 = 0.2 or 20%)
// - Today: 1-year USD/BRL implied volatility is harder to get exact without options data, but 2025 estimates suggest 15-20%. We’ll lower it to reflect a calmer market.
// - New Value: v0 = 0.0225 (volatility = square root of 0.0225 = 0.15 or 15%).

// 7. kappa - Mean Reversion Rate
// - Original: 2.0 (variance reverts in ~0.5 years)
// - This is a model parameter, not tied to today’s spot data. It’s reasonable for FX, so we’ll keep it.
// - New Value: kappa = 2.0 (no change).

// 8. theta - Long-Term Variance
// - Original: 0.04 (long-term volatility = square root of 0.04 = 0.2 or 20%)
// - Should match long-term expectations. USD/BRL volatility historically is 15-25%, so we’ll align with v0’s lower estimate.
// - New Value: theta = 0.0225 (long-term volatility = 15%).

// 9. sigma - Volatility of Variance
// - Original: 0.3 (30% vol-of-vol)
// - Controls variance fluctuations. Without options data to calibrate, 0.3 is typical for FX. We’ll keep it.
// - New Value: sigma = 0.3 (no change).

// 10. rho - Correlation
// - Original: -0.7 (negative correlation between spot and volatility)
// - Common in FX (volatility rises as currency weakens). No reason to change without skew data.
// - New Value: rho = -0.7 (no change).


@Service
public class HestonModelTrapezoidalIntegration {
    //TrapezoidalIntegration
    // Market and option parameters
    private final double spotPrice;          // S0 - Spot Price (USD/BRL)
    private final double strikePrice;         // K  - Strike Price
    private final double timeToMaturity;     // T  - Time to Maturity (years)
    private final double usdRiskFreeRate;    // r  - USD Risk-Free Rate
    private final double brlRiskFreeRate;    // q  - BRL Risk-Free Rate
    private final double initialVariance;    // v0 - Initial Variance
    private final double meanReversionRate; // kappa - Mean Reversion Rate
    private final double longTermVariance;   // theta - Long-Term Variance
    private final double volOfVariance;      // sigma - Volatility of Variance
    private final double correlation;        // rho   - Correlation

    /**
     * Constructor for the Heston Model.
     *
     * @param spotPrice          Spot price of the underlying asset.
     * @param strikePrice         Strike price of the option.
     * @param timeToMaturity     Time to maturity of the option in years.
     * @param usdRiskFreeRate    Risk-free interest rate in USD.
     * @param brlRiskFreeRate    Risk-free interest rate in BRL.
     * @param initialVariance    Initial variance of the asset.
     * @param meanReversionRate Mean reversion rate of the variance.
     * @param longTermVariance   Long-term variance.
     * @param volOfVariance      Volatility of the variance.
     * @param correlation        Correlation between the asset and its variance.
     */
    public HestonModelTrapezoidalIntegration(double spotPrice,
                                             double strikePrice,
                                             double timeToMaturity,
                                             double usdRiskFreeRate,
                                             double brlRiskFreeRate,
                                             double initialVariance,
                                             double meanReversionRate,
                                             double longTermVariance,
                                             double volOfVariance,
                                             double correlation) {
        this.spotPrice = spotPrice;
        this.strikePrice = strikePrice;
        this.timeToMaturity = timeToMaturity;
        this.usdRiskFreeRate = usdRiskFreeRate;
        this.brlRiskFreeRate = brlRiskFreeRate;
        this.initialVariance = initialVariance;
        this.meanReversionRate = meanReversionRate;
        this.longTermVariance = longTermVariance;
        this.volOfVariance = volOfVariance;
        this.correlation = correlation;
    }

    // Complex number class for characteristic function
    static class Complex {
        double real, imag;

        Complex(double real, double imag) {
            this.real = real;
            this.imag = imag;
        }

        Complex add(Complex other) {
            return new Complex(this.real + other.real, this.imag + other.imag);
        }

        Complex subtract(Complex other) {
            return new Complex(this.real - other.real, this.imag - other.imag);
        }

        Complex multiply(Complex other) {
            return new Complex(
                this.real * other.real - this.imag * other.imag,
                this.real * other.imag + this.imag * other.real
            );
        }

        Complex scale(double scalar) {
            return new Complex(this.real * scalar, this.imag * scalar);
        }

        Complex exp() {
            double expReal = Math.exp(this.real);
            return new Complex(expReal * Math.cos(this.imag), expReal * Math.sin(this.imag));
        }

        double real() {
            return this.real;
        }

        public Complex log() {
            double r = Math.sqrt(this.real * this.real + this.imag * this.imag);
            double theta = Math.atan2(this.imag, this.real);
            return new Complex(Math.log(r), theta);
        }

        public Complex divide(Complex other) {
            double denom = other.real * other.real + other.imag * other.imag;
            double safeDenom = Math.max(denom, 1e-10);
            double real = (this.real * other.real + this.imag * other.imag) / safeDenom;
            double imag = (this.imag * other.real - this.real * other.imag) / safeDenom;
            return new Complex(real, imag);
        }

        // Placeholder for complex square root (simplified)
        public Complex sqrtComplex() {
            double r = Math.sqrt(Math.sqrt(this.real * this.real + this.imag * this.imag));
            double theta = Math.atan2(this.imag, this.real) / 2;
            if (Double.isNaN(r) || Double.isNaN(theta)) {
                return new Complex(0, 0); // Fallback
            }
            return new Complex(r * Math.cos(theta), r * Math.sin(theta));
        }


    }

    // Heston characteristic function
    private Complex characteristicFunction(double phi, int j) {
        Complex i = new Complex(0, 1);
        Complex u = (j == 1) ? new Complex(0.5, 0) : new Complex(-0.5, 0);
        double a = meanReversionRate * longTermVariance;
        Complex b = (j == 1) ? new Complex(meanReversionRate - correlation * volOfVariance, 0) : new Complex(meanReversionRate, 0);
        Complex phiI = i.scale(phi);

        Complex bMinusRhoSigmaPhiI = b.subtract(phiI.scale(correlation * volOfVariance));
        Complex dTerm = bMinusRhoSigmaPhiI.multiply(bMinusRhoSigmaPhiI)
                .subtract(new Complex(-phi * phi, 0).add(u.multiply(phiI).scale(2))
                        .scale(volOfVariance * volOfVariance));
        Complex d = dTerm.sqrtComplex();

        Complex gNumer = bMinusRhoSigmaPhiI.add(d);
        Complex gDenom = bMinusRhoSigmaPhiI.subtract(d);
        double gDenomReal = gDenom.real();
        double safeGDenom = (Math.abs(gDenomReal) < 1e-10) ? ((gDenomReal < 0) ? -1e-10 : 1e-10) : gDenomReal;
        Complex g = gNumer.scale(1.0 / safeGDenom);

        Complex dT = d.scale(timeToMaturity);
        Complex expDT = new Complex(Math.cos(dT.imag), Math.sin(dT.imag)).scale(Math.exp(dT.real()));
        Complex one = new Complex(1, 0);
        Complex geDT = g.multiply(expDT);
        Complex term1 = one.subtract(geDT);
        Complex term2 = one.subtract(g);
        Complex logArg = term1.divide(term2); // Use Complex divide
        Complex lnTerm = logArg.log();

        Complex C = phiI.scale((usdRiskFreeRate - brlRiskFreeRate) * timeToMaturity)
                .add(bMinusRhoSigmaPhiI.add(d).scale(timeToMaturity).subtract(lnTerm.scale(2))
                        .scale(a / (volOfVariance * volOfVariance)));

        Complex Dnumer = bMinusRhoSigmaPhiI.add(d).multiply(one.subtract(expDT));
        Complex Ddenom = term1.scale(volOfVariance * volOfVariance);
        Complex D = Dnumer.divide(Ddenom);

        Complex exponent = C.add(D.scale(initialVariance)).add(phiI.scale(Math.log(spotPrice)));
        double expReal = exponent.real();
        if (expReal > 700) expReal = 700; // Prevent overflow
        if (expReal < -700) expReal = -700; // Prevent underflow
        double expMag = Math.exp(expReal);
        return new Complex(expMag * Math.cos(exponent.imag), expMag * Math.sin(exponent.imag));
    }



    private void testCharacteristicFunction() {
        System.out.println("BEGIN testCharacteristicFunction");
        double[] phiValues = {0, 1, 10};
        for (double phi : phiValues) {
            for (int j = 1; j <= 2; j++) {
                characteristicFunction(phi, j);
            }
        }
        System.out.println("END testCharacteristicFunction");
    }

    private double integrate(double phiMax, int nSteps, int j) {
        double deltaPhi = phiMax / nSteps;
        double sum = 0.0;
        for (int i = 0; i <= nSteps; i++) {
            double phi = i * deltaPhi;
            Complex cf = characteristicFunction(phi, j);
            double integrand;
            if (phi == 0.0) {
                integrand = 0.0;
            } else {
                Complex expTerm = new Complex(Math.cos(-phi * Math.log(strikePrice)), Math.sin(-phi * Math.log(strikePrice)));
                Complex numerator = expTerm.multiply(cf).multiply(new Complex(0, -1));
                integrand = numerator.real() / phi; // Revert sign flip
                System.out.println("j: " + j + ", phi: " + phi + ", cf: " + cf.real + " + " + cf.imag + "i, integrand: " + integrand);
            }
            if (Double.isNaN(integrand) || Double.isInfinite(integrand)) {
                integrand = 0.0;
            }
            sum += (i == 0 || i == nSteps) ? integrand / 2 : integrand;
        }
        double result = sum * deltaPhi / Math.PI;
        System.out.println("P" + j + " integral contribution: " + result);
        return result;
    }

    // Calculate call option price
    public double callPrice() {
        double phiMax = 200; // Wider range
        int nSteps = 20000;  // More steps
        double P1 = 0.5 + integrate(phiMax, nSteps, 1);
        double P2 = 0.5 + integrate(phiMax, nSteps, 2);
        System.out.println("P1: " + P1 + ", P2: " + P2);
        double price = spotPrice * Math.exp(-brlRiskFreeRate * timeToMaturity) * P1 - strikePrice * Math.exp(-usdRiskFreeRate * timeToMaturity) * P2;
        return price > 0 ? price : 0; // Ensure non-negative price
    }

    private static double getCallPriceValue() {
        /*
Original Parameters:
S0=5.80, K=5.85, T=1.0, r=0.045, q=0.105, v0=0.04, kappa=2.0, theta=0.04, sigma=0.3, rho=-0.7
Price = 0.1998

Updated Parameters (Feb 21, 2025):
S0=5.75, K=5.85, T=1.0, r=0.0425, q=0.1075, v0=0.0225, kappa=2.0, theta=0.0225, sigma=0.3, rho=-0.7
Expected Price = ~0.12 (lower due to 15% volatility vs. 20%)

        * */
        double S0 = 5.80;    // Reverted spot price
        double K = 5.85;     // Strike price (unchanged)
        double T = 1.0;      // Time to maturity (unchanged)
        double r = 0.0425;   // Updated USD risk-free rate
        double q = 0.1075;   // Updated BRL risk-free rate
        double v0 = 0.0225;  // Updated initial variance
        double kappa = 2.0;  // Unchanged
        double theta = 0.0225; // Updated long-term variance
        double sigma = 0.3;  // Unchanged
        double rho = -0.7;   // Unchanged

        HestonModelTrapezoidalIntegration heston = new HestonModelTrapezoidalIntegration(S0, K, T, r, q, v0, kappa, theta, sigma, rho);
        //heston.testCharacteristicFunction(); // Debug samples
        return heston.callPrice();
    }

    public static void main(String[] args) {
        double callPrice = getCallPriceValue();
        System.out.printf("Heston Call Option Price for BRL/USD: %.4f%n", callPrice);

    }
}
/*

Parameters

S0 = 5.80 (spot price, BRL/USD)
K = 5.85 (strike price)
T = 1.0 (time to maturity, years)
r = 0.045 (USD risk-free rate)
q = 0.105 (BRL risk-free rate)
v = sqrt(v0) = sqrt(0.04) = 0.2 (volatility, using initial variance from Heston)

Compare with Black-Scholes Formula

    C = S0 * e^(-q * T) * N(d1) - K * e^(-r * T) * N(d2)
    d1 = [ln(S0 / K) + (r - q + 0.5 * v^2) * T] / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    N(x) = cumulative normal distribution function (approximated below)

1. Compute d1

    ln(S0 / K) = ln(5.80 / 5.85) = ln(0.991453) ~= -0.008584
    r - q = 0.045 - 0.105 = -0.06
    0.5 * v^2 = 0.5 * 0.2^2 = 0.5 * 0.04 = 0.02
    (r - q + 0.5 * v^2) * T = (-0.06 + 0.02) * 1 = -0.04
    Numerator = ln(S0 / K) + (r - q + 0.5 * v^2) * T = -0.008584 + (-0.04) = -0.048584
    Denominator = v * sqrt(T) = 0.2 * sqrt(1) = 0.2
    d1 = -0.048584 / 0.2 = -0.24292

2. Compute d2
    d2 = d1 - v * sqrt(T) = -0.24292 - 0.2 * 1 = -0.44292

3. Approximate N(d1) and N(d2)

    N(x) ~= 1 / (1 + e^(-1.5955 * x))  (sigmoid-based, rough but ASCII-friendly)

    More precise methods exist (e.g., polynomial approximations), but this keeps it simple.

    For d1 = -0.24292:

    -1.5955 * (-0.24292) = 0.38758
    e^(0.38758) ~= 1.4735
    1 + e^(0.38758) = 1 + 1.4735 = 2.4735
    N(d1) = 1 / 2.4735 ~= 0.4042

    For d2 = -0.44292:

    -1.5955 * (-0.44292) = 0.70669
    e^(0.70669) ~= 2.0272
    1 + e^(0.70669) = 1 + 2.0272 = 3.0272
    N(d2) = 1 / 3.0272 ~= 0.3303

4. Compute Discount Factors

    e^(-q * T) = e^(-0.105 * 1) = e^(-0.105) ~= 0.9003
    e^(-r * T) = e^(-0.045 * 1) = e^(-0.045) ~= 0.9560

5. Compute Call Price

    C = S0 * e^(-q * T) * N(d1) - K * e^(-r * T) * N(d2)
    First term = 5.80 * 0.9003 * 0.4042 ~= 2.1098
    Second term = 5.85 * 0.9560 * 0.3303 ~= 1.8478
    C = 2.1098 - 1.8478 ~= 0.2620

    Result
        Black-Scholes Call Price: ~0.2620

6. Summary Black-Scholes:
    d1 = -0.24292, d2 = -0.44292
    N(d1) ~= 0.4042, N(d2) ~= 0.3303
    C = 5.80 * 0.9003 * 0.4042 - 5.85 * 0.9560 * 0.3303 ~= 0.2620

    Heston (your result): 0.1998

Notes
1. Approximation: The N(x) approximation here is crude. Using a standard normal CDF table or better formula (e.g., Abramowitz-Stegun):
    * 𝑁(−0.24292)≈0.404≈0.404 (close to our 0.4042)
    * 𝑁(−0.44292)≈0.329≈0.329 (close to our 0.3303)
    This confirms the rough calculation is reasonable, but a precise N(x) might adjust the price slightly (e.g., closer to 0.26).

2. Why Heston Differs:
    *Black-Scholes assumes constant volatility (20% here), while Heston models stochastic volatility (𝑣0=0.04, 𝜃=0.04, 𝜎=0.3, ρ=−0.7).

    *The negative correlation (ρ=−0.7) in Heston reduces the option price for out-of-the-money calls by skewing volatility downward when the spot falls, unlike Black-Scholes’ flat volatility.

    *Your Heston price (0.1998) being lower than Black-Scholes (0.2620) aligns with this effect.

3. Correction Note: In my prior Black-Scholes estimate (~0.155), I misstepped—rechecking with precise CDF values gives ~0.26, consistent here. Apologies for that!

*/


