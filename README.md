# DoubleDoubleStatistic
 Double-Double Statistic Implements

## Requirement
.NET 8.0  
[DoubleDouble](https://github.com/tk-yoshimura/DoubleDouble)  
[DoubleDoubleComplex](https://github.com/tk-yoshimura/DoubleDoubleComplex)  

## Install
[Download DLL](https://github.com/tk-yoshimura/DoubleDoubleStatistic/releases)  
[Download Nuget](https://www.nuget.org/packages/tyoshimura.DoubleDouble.Statistic/)  

## Implemented Distributions

### Continuous

| category   | distribution      | PDF      | CDF      | quantile | statistic | fitting  | random generation | note                                             |
| ---------- | ----------------- | -------- | -------- | -------- | --------- | -------- | ----------------- | ------------------------------------------------ |
| stable     | cauchy            | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | delta             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | holtsmark         | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | landau            | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | levy              | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | map-airy          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | normal            | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | sas point5        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
| linearity  | cosine            | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | davis             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | CDF and Quantile take longer to calculate.       |
|            | frechet           | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | gumbel            | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | johnson sb        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | johnson su        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | laplace           | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | logistic          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | skew normal       | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | Quantile take longer to calculate.               |
|            | uniform           | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | u quadratic       | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | weibull           | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
| scalable   | birnbaum saunders | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | exponential       | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | folded normal     | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | gamma             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | gompertz          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | half normal       | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | Quantile take longer to calculate.               |
|            | hyperbolic secant | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | inverse gauss     | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | Quantile take longer to calculate.               |
|            | log logistic      | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | lomax             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | maxwell           | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | pareto            | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | rayleigh          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | voigt             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | CDF and Quantile take longer to calculate.       |
|            | wigner semicircle | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
| continuous | alpha             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | arcsine           | &#x2714; | &#x2714; | &#x2714; | &#x2714;  |          | &#x2714;          |                                                  |
|            | argus             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          |                                                  |
|            | benktander        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | Quantile take longer to calculate.               |
|            | beta              | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | beta prime        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | bradford          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | burr              | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | chi               | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | chi square        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | dagum             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | fisher z          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | fisk              | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | inverse gamma     | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | inverse chi       | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | inverse chi sq    | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | irwin hall        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  |          | &#x2714;          | n &leq; 128                                      |
|            | kumaraswamy       | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | log normal        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | nakagami          | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | noncentral chi sq | &#x2714; | &#x2714; | &#x2714; | &#x2714;  |          | &#x2714;          | Accuracy decreases when non-centricity is large. |
|            | noncentral f      | &#x2714; | &#x2714; | &#x2714; | &#x2714;  |          | &#x2714;          | Accuracy decreases when non-centricity is large. |
|            | noncentral t      | &#x2714; | &#x2714; | &#x2714; | &#x2714;  |          | &#x2714;          | Accuracy decreases when non-centricity is large. |
|            | power             | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | reciprocal        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | rice              | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x26A0; | &#x26A0;          | CDF and Quantile take longer to calculate.       |
|            | snedecor f        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | student t         | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |
|            | triangular        | &#x2714; | &#x2714; | &#x2714; | &#x2714;  | &#x2714; | &#x2714;          |                                                  |

### Discrete

| category | distribution      | PMF      | statistic | note |
| -------- | ----------------- | -------- | --------- | ---- |
| discrete | bernoulli         | &#x2714; | &#x2714;  |      |
|          | binary            | &#x2714; | &#x2714;  |      |
|          | binomial          | &#x2714; | &#x2714;  |      |
|          | categorical       | &#x2714; | &#x2714;  |      |
|          | discrete uniform  | &#x2714; | &#x2714;  |      |
|          | gausskuzmin       | &#x2714; | &#x2714;  |      |
|          | geometric         | &#x2714; | &#x2714;  |      |
|          | hyper geometric   | &#x2714; | &#x2714;  |      |
|          | logarithmic       | &#x2714; | &#x2714;  |      |
|          | negative binomial | &#x2714; | &#x2714;  |      |
|          | poisson           | &#x2714; | &#x2714;  |      |
|          | skellam           | &#x2714; | &#x2714;  |      |
|          | yule simon        | &#x2714; | &#x2714;  |      |
|          | zipf              | &#x2714; | &#x2714;  |      |

## Usage

```cs
NormalDistribution dist = new(mu: 1, sigma: 3);

// PDF
for (ddouble x = -4; x <= 4; x += 0.125) {
    ddouble pdf = dist.PDF(x);

    Console.WriteLine($"pdf({x})={pdf}");
}

// CDF
for (ddouble x = -4; x <= 4; x += 0.125) {
    ddouble ccdf = dist.CDF(x, Interval.Upper);

    Console.WriteLine($"ccdf({x})={ccdf}");
}

// Quantile
for (int i = 0; i <= 10; i++) {
    ddouble p = (ddouble)i / 10;
    ddouble x = dist.Quantile(p, Interval.Upper);

    Console.WriteLine($"cquantile({p})={x}");
}

// Statistic
Console.WriteLine($"Support={dist.Support}");
Console.WriteLine($"Mu={dist.Mu}");
Console.WriteLine($"Sigma={dist.Sigma}");
Console.WriteLine($"Mean={dist.Mean}");
Console.WriteLine($"Median={dist.Median}");
Console.WriteLine($"Mode={dist.Mode}");
Console.WriteLine($"Variance={dist.Variance}");
Console.WriteLine($"Skewness={dist.Skewness}");
Console.WriteLine($"Kurtosis={dist.Kurtosis}");
Console.WriteLine($"Entropy={dist.Entropy}");

// Random Sampling
Random random = new(1234);
double[] xs = dist.Sample(random, 100000).ToArray();
```

## Typical parameter symbols
| category                 | symbol       | note                        |
| ------------------------ | ------------ | --------------------------- |
| support parameter        | k            |                             |
|                          | a, b         | uniform                     |
|                          | a, b, c      | triangular                  |
| shape parameter          | alpha        |                             |
|                          | alpha, beta  | beta, beta prime            |
|                          | gamma, delta | johnson sb, su              |
|                          | eta          | gompertz                    |
|                          | nu           | chi, chisq, student t       |
|                          | n            | irwin hall                  |
|                          | n, m         | fisher z, snedecor f        |
|                          | c            | stable distributions        |
| location parameter       | mu           |                             |
| scale parameter          | sigma        | error-related distributions |
|                          | theta        | time-related distributions  |
|                          | s, r         | otherwise                   |
| non-centricity parameter | lambda       |                             |
|                          | mu           | non-central student t       |

## To Be Updated

### Not Implemented Distributions

(N/A)

### Not Implemented Functions

- Fitting (Discrete)
- Statistic Test

## Licence
[MIT](https://github.com/tk-yoshimura/DoubleDoubleStatistic/blob/main/LICENSE)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)
