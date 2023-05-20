package com.nnonaka.option.price.controller;

import com.nnonaka.option.price.core.model.EntryProtocol;
import com.nnonaka.option.price.core.model.OutProtocol;
import com.nnonaka.option.price.service.GarmanKohlhagenOptionFFTPricing;
import com.nnonaka.option.price.service.GarmanKohlhagenOptionPricing;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class GarmanKohlhagenController {
    @Autowired
    GarmanKohlhagenOptionPricing garmanKohlhagenOptionPricing;
    @Autowired
    GarmanKohlhagenOptionFFTPricing garmanKohlhagenOptionFFTPricing;

    @RequestMapping(value = "/v1/option/garman/simulation", method = RequestMethod.POST, consumes = MediaType.APPLICATION_JSON_VALUE)
    public OutProtocol priceSimulation(@RequestBody final EntryProtocol entryProtocol){
        double spotPrice = 5.25; // Spot price of the currency pair (USD/BRL)
        double strikePrice = 5.20; // Strike price
        double riskFreeRate = 0.05; // Risk-free interest rate
        double volatility = 0.15; // Volatility of the currency pair
        double timeToExpiration = 1.0; // Time to expiration in years

        String callPriceResult = garmanKohlhagenOptionPricing.europeanCallOptionPrice(spotPrice, strikePrice, riskFreeRate, volatility, timeToExpiration);
        OutProtocol outProtocol= new OutProtocol();
        String data= String.format("%s", callPriceResult);
        outProtocol.setData(data);
        return outProtocol;

    }
    @RequestMapping(value = "/v1/option/garman/usdbrl", method = RequestMethod.POST, consumes = MediaType.APPLICATION_JSON_VALUE)
    public OutProtocol priceUsdBrlSimulation(@RequestBody final EntryProtocol entryProtocol){
        double spotPrice = 5.25; // Spot price of the currency pair (USD/BRL)
        double strikePrice = 5.20; // Strike price
        double brazilRiskFreeRate = 0.116; // Risk-free interest rate in Brazil
        double usaRiskFreeRate = 0.0525; // Risk-free interest rate in the USA
        double volatility = 0.15; // Volatility of the currency pair
        double timeToExpiration = 0.083; // Time to expiration in years

        String callPrice = garmanKohlhagenOptionPricing.europeanCallOptionPrice(spotPrice, strikePrice, brazilRiskFreeRate, usaRiskFreeRate, volatility, timeToExpiration);
        OutProtocol outProtocol= new OutProtocol();
        String data= String.format("%s", callPrice);
        outProtocol.setData(data);
        return outProtocol;

    }
    @RequestMapping(value = "/v1/option/garman/fft/simulation", method = RequestMethod.POST, consumes = MediaType.APPLICATION_JSON_VALUE)
    public OutProtocol priceUsdBrlFftSimulation(@RequestBody final EntryProtocol entryProtocol){
        double spotPrice = 5.25; // Spot price of the currency pair (USD/BRL)
        double strikePrice = 5.20; // Strike price
        double riskFreeRate = 0.05; // Risk-free interest rate
        double volatility = 0.15; // Volatility of the currency pair
        double timeToExpiration = 0.083; // Time to expiration in years

        String callPrice = garmanKohlhagenOptionFFTPricing.europeanCallOptionPrice(spotPrice, strikePrice, riskFreeRate, volatility, timeToExpiration);
        OutProtocol outProtocol= new OutProtocol();
        String data= String.format("%s", callPrice);
        outProtocol.setData(data);
        return outProtocol;

    }

}
