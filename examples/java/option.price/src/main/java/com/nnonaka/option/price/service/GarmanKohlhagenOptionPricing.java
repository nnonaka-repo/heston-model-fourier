package com.nnonaka.option.price.service;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.springframework.stereotype.Service;

@Service
public class GarmanKohlhagenOptionPricing {
    public String europeanCallOptionPrice(double spotPrice, double strikePrice, double riskFreeRate, double volatility, double timeToExpiration) {

        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + 0.5 * Math.pow(volatility, 2)) * timeToExpiration) / (volatility * Math.sqrt(timeToExpiration));
        double d2 = d1 - volatility * Math.sqrt(timeToExpiration);

        NormalDistribution distribution = new NormalDistribution();
        double N1 = distribution.cumulativeProbability(d1);
        double N2 = distribution.cumulativeProbability(d2);

        double callPrice = spotPrice * Math.exp(-riskFreeRate * timeToExpiration) * N1 - strikePrice * Math.exp(-riskFreeRate * timeToExpiration) * N2;
        return  String.format("callPrice : %f", callPrice);
    }
    public String europeanCallOptionPrice(double spotPrice, double strikePrice, double brazilRiskFreeRate, double usaRiskFreeRate, double volatility, double timeToExpiration) {
        double d1 = (Math.log(spotPrice / strikePrice) + (brazilRiskFreeRate - usaRiskFreeRate + 0.5 * Math.pow(volatility, 2)) * timeToExpiration) / (volatility * Math.sqrt(timeToExpiration));
        double d2 = d1 - volatility * Math.sqrt(timeToExpiration);

        NormalDistribution distribution = new NormalDistribution();
        double N1 = distribution.cumulativeProbability(d1);
        double N2 = distribution.cumulativeProbability(d2);

        double callPrice = spotPrice * Math.exp(-usaRiskFreeRate * timeToExpiration) * N1 - strikePrice * Math.exp(-brazilRiskFreeRate * timeToExpiration) * N2;
        return  String.format("callPrice : %f", callPrice);
    }
}
