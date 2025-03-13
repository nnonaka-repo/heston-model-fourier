package com.nnonaka.option.price.service;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class GarmanKohlhagenOptionFFTOptionPricingClaude35 {

    private static final int N = 4096; // Número de pontos da FFT (deve ser uma potência de 2)
    private static final double alpha = 1.5; // Parâmetro de amortecimento

    public static double[] priceOptions(double S, double r, double q, double sigma, double T, double[] strikes, boolean isPut) {
        double eta = 0.05; // Espaçamento da grade. Espaços entre os pontos
        double lambda = 2 * Math.PI / (N * eta);
        double b = lambda * N / 2;

        Complex[] x = new Complex[N];
        for (int k = 0; k < N; k++) {
            double v = k * eta;
            Complex u = new Complex(v, -(alpha + 1));
            Complex psi = characteristicFunction(u, S, r, q, sigma, T);
            Complex denominator = u.multiply(u).add(Complex.I.multiply(u)).subtract(Math.pow(alpha,2) + alpha);
            Complex modifiedPsi = psi.divide(denominator);
            x[k] = modifiedPsi.multiply(Math.exp(-r * T)).multiply(Complex.I.negate());
        }

        FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] y = transformer.transform(x, TransformType.FORWARD);

        double[] callPrices = new double[strikes.length];
        double[] putPrices = new double[strikes.length];
        for (int i = 0; i < strikes.length; i++) {
            double k = Math.log(strikes[i]);
            int index = (int) Math.round((k + b) / lambda);
            if (index >= 0 && index < N) {
                double simpson = (index == 0 || index == N - 1) ? 1 : (2 + 2 * (index % 2)) / 3.0;
                callPrices[i] = Math.exp(-alpha * k) * y[index].getReal() / Math.PI * simpson * eta;
                putPrices[i] = callPrices[i] - S * Math.exp(-q * T) + strikes[i] * Math.exp(-r * T);
            }
        }

        return isPut?putPrices:callPrices;
    }

    private static Complex characteristicFunction(Complex u, double S, double r, double q, double sigma, double T) {
        Complex i = Complex.I;
        //Complex mu = new Complex(riskFreeRate - dividend - 0.5 * Math.pow(volatility, 2));
        Complex mu = new Complex(r - q - 0.5 * Math.pow(sigma,2));//drift
        Complex exponent = i.multiply(u).multiply(Math.log(S) + mu.getReal() * T)
                .subtract(u.multiply(u).multiply(0.5 * Math.pow(sigma,2) * T));
        return exponent.exp();
    }

/*
double q = 0.02; // Taxa de dividendos de 2% ao ano

Para ações: q = Dividendo anual total / Preço atual da ação
Para índices: Pode ser calculado como a média ponderada dos rendimentos de dividendos das ações componentes.

Para ações que não pagam dividendos, q = 0.
Para opções sobre moedas, q pode representar a taxa de juros estrangeira.
Em commodities, pode representar o custo de armazenamento ou o convenience yield.
*/
/*
Eta (η):

No contexto da FFT para precificação de opções, eta é um parâmetro de discretização.

Representa o espaçamento entre os pontos na grade de integração numérica.

É usado para determinar a resolução da FFT e, consequentemente, a precisão do cálculo dos preços das opções.
 * */
/*
* Definição Delta bump:

Delta bump é uma técnica de análise de sensibilidade usada para calcular o delta de uma opção ou outro derivativo financeiro.
Envolve "perturbar" ou "bumpar" (alterar ligeiramente) o preço do ativo subjacente e observar como isso afeta o preço do derivativo.
Processo básico:

Calcular o preço do derivativo com o preço atual do ativo subjacente.
Recalcular o preço do derivativo após um pequeno aumento no preço do ativo subjacente.
A diferença entre esses dois preços, dividida pela mudança no preço do ativo, fornece uma aproximação do delta.
* */
public static void main(String[] args) {
    double S = 5057.8;          // Taxa de câmbio spot USD/BRL (1 dólar = 5.6 reais)
    double rDomestica = 0.045;    // Taxa de juros americana (anualizada)
    double rEstrangeira = 0.10;   // Taxa de juros brasileira (anualizada)
    double sigma = 0.10;      // Volatilidade de 10% (anualizada)
    double T = 1.0 / 12;      // Tempo até o vencimento: 1 mês (como fração do ano)

    double[] strikes = {5019.867, 5057.8, 5095.734, 5800, 6000}; // Preços de exercício

        double[] callPrices = priceOptions(S, rDomestica, rEstrangeira, sigma, T, strikes,false);
        System.out.println("Preços das opções de compra usando FFT:");
        for (int i = 0; i < strikes.length; i++) {
            System.out.printf("Strike: %.2f, Preço: %.4f%n", strikes[i], callPrices[i]);
        }

        double[] putPrices = priceOptions(S, rDomestica, rEstrangeira, sigma, T, strikes,true);
        System.out.println("\n\rPreços das opções de venda (PUT) usando FFT:");
        for (int i = 0; i < strikes.length; i++) {
            System.out.printf("Strike: %.2f, Preço: %.4f%n", strikes[i], putPrices[i]);
        }
    }
}

/*

Exemplo lúdico sobre esse drift

Complex mu = new Complex(riskFreeRate - dividend - 0.5 * Math.pow(volatility, 2));

* Imagine um parque aquático chamado "Mercado Financeiro" com um rio artificial chamado "Rio dos Ativos". Neste rio, as pessoas flutuam em boias coloridas, cada uma representando um ativo diferente (ações, moedas, commodities, etc.).

O Rio e suas Correntes:

O rio tem uma corrente principal, que representa o drift.
A direção e a força da corrente são o equivalente à tendência geral do mercado.
As Boias (Ativos):

Cada boia (ativo) está sujeita à corrente principal (drift) do rio.
Algumas boias são mais pesadas (ações de empresas estáveis), outras mais leves (ações voláteis).
O Drift como Corrente:

Se o drift é positivo, a corrente do rio flui suavemente para frente (mercado em alta).
Se o drift é negativo, algumas partes do rio têm contracorrentes (mercado em baixa).
Um drift zero seria como um lago calmo, sem corrente perceptível.
Elementos que Afetam o Drift:

"Bombas de Água" (taxa de juros): Controlam a força da corrente principal.
"Cachoeiras" (pagamento de dividendos): Criam quedas no nível da água, afetando a velocidade das boias.
Volatilidade como Ondas:

Além da corrente principal, há ondas no rio (volatilidade).
Algumas partes do rio são mais agitadas, outras mais calmas.
Precificação de Opções com FFT:

Imagine um fotógrafo mágico (o algoritmo FFT) que pode tirar uma foto do rio inteiro de uma só vez.
Esta foto não só mostra onde as boias estão agora, mas prevê onde elas provavelmente estarão no futuro.
O Papel do Drift na FFT:

O drift (corrente do rio) ajuda o fotógrafo a prever para onde as boias estão se movendo em geral.
Sem considerar o drift, a foto ficaria borrada ou imprecisa.
Ajuste de Risco:

Há um "modo seguro" no parque onde todas as boias se movem com a mesma velocidade média (precificação neutra ao risco).
Neste modo, a corrente é ajustada para (r - q), onde r é a velocidade padrão do parque e q é o "atrito" causado pelas cachoeiras.
Exemplo Prático:

Se o drift é positivo (corrente para frente), as boias tendem a terminar mais à frente no rio.
Isso é como opções de compra (calls) se tornando mais valiosas em um mercado em alta.
Desafio do Fotógrafo (FFT):

O fotógrafo precisa considerar a corrente (drift), as ondas (volatilidade), e os obstáculos (dividendos, taxas de juros) para tirar uma foto precisa do futuro do rio.
* */