package com.nnonaka.option.price.controller;

import com.nnonaka.option.price.core.model.EntryProtocol;
import com.nnonaka.option.price.core.model.OutProtocol;
import com.nnonaka.option.price.service.CarrMadanOptionPricingFFTService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class CarrMadanController {

    @Autowired
    CarrMadanOptionPricingFFTService carrMadanOptionPricingService;
    @RequestMapping(value = "/v1/option/simulation", method = RequestMethod.POST, consumes = MediaType.APPLICATION_JSON_VALUE)
    public OutProtocol priceSimulation(@RequestBody final EntryProtocol entryProtocol){
        double S = 100.0; // underlying asset price
        double K = 100.0; // strike price
        double r = 0.05; // risk-free interest rate
        double sigma = 0.2; // volatility
        double T = 1.0; // time to expiration
        String callPriceResult = carrMadanOptionPricingService.europeanCallOptionPrice(K, sigma,S, r, T);
        OutProtocol outProtocol= new OutProtocol();
        String data= String.format("%s", callPriceResult);
        outProtocol.setData(data);
        return outProtocol;
    }
}
