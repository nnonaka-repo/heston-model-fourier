package com.nnonaka.option.price.service;

public class GarmanKohlhagenSemFFT {
    private static double normalCDF(double x) {
        double k = 1.0 / (1.0 + 0.2316419 * Math.abs(x));
        double y = ((((1.330274429 * k - 1.821255978) * k + 1.781477937) *
                k - 0.356563782) * k + 0.319381530) * k;
        y = 1.0 - 1.0 / Math.sqrt(2 * Math.PI) * Math.exp(-x * x / 2.0) * y;
        return x < 0 ? 1.0 - y : y;
    }

    // Função para calcular d1 e d2
    private static double[] calculateD1D2(double S, double K, double r_d, double r_f, double sigma, double T) {
        double d1 = (Math.log(S / K) + (r_d - r_f + 0.5 * sigma * sigma) * T) / (sigma * Math.sqrt(T));
        double d2 = d1 - sigma * Math.sqrt(T);
        return new double[]{d1, d2};
    }

    // Função para calcular o preço da opção de compra (call)
    public static double callPrice(double S, double K, double r_d, double r_f, double sigma, double T) {
        double[] d = calculateD1D2(S, K, r_d, r_f, sigma, T);
        double d1 = d[0];
        double d2 = d[1];

        return S * Math.exp(-r_f * T) * normalCDF(d1) - K * Math.exp(-r_d * T) * normalCDF(d2);
    }

    // Função para calcular o preço da opção de venda (put)
    public static double putPrice(double S, double K, double r_d, double r_f, double sigma, double T) {
        double[] d = calculateD1D2(S, K, r_d, r_f, sigma, T);
        double d1 = d[0];
        double d2 = d[1];

        return K * Math.exp(-r_d * T) * normalCDF(-d2) - S * Math.exp(-r_f * T) * normalCDF(-d1);
    }

    public static void main(String[] args) {
        // Exemplo de uso
        double S = 5.5000;  // Taxa de câmbio spot EUA/BRL
        double K = 5.6500;  // Preço de exercício
        double r_d = 0.04;  // Taxa de juros doméstica (USD)
        double r_f = 10.05;  // Taxa de juros estrangeira (BRL)
        double sigma = 10.15; // Volatilidade
        double T = 0.083;     // Tempo até o vencimento 1 mês (em anos)

        double callOption = callPrice(S, K, r_d, r_f, sigma, T);
        double putOption = putPrice(S, K, r_d, r_f, sigma, T);

        System.out.printf("Preço da opção de compra (Call): %.4f%n", callOption);
        System.out.printf("Preço da opção de venda (Put): %.4f%n", putOption);
    }
}
