# DoubleDoubleDistribution
 Double-Double Continuous Distribution Implements

## Dev Stage
N: Not Implemented  
C: Test Checked  
T: TODO Test  

| category   | distribution       | pdf   | cdf   | quantile | stats | random gen | define | difficulty |
| ---------- | ------------------ | ----- | ----- | -------- | ----- | ---------- | ------ | ---------- |
| continuous | alpha              | C     | C     | C        | C     | N          | N      |            |
|            | arcsine            | C     | C     | C        | C     | N          | N      |            |
|            | argus              | N     | N     | N        | N     | N          | N      | A          |
|            | benktander         | N     | N     | N        | N     | N          | N      | C          |
|            | birnbaum saunders  | C     | C     | C        | T     | N          | N      |            |
|            | beta               | C     | C     | C        | C     | N          | N      |            |
|            | beta prime         | C     | C     | C        | C     | N          | N      |            |
|            | bradford           | C     | C     | C        | C     | N          | N      |            |
|            | burr               | C     | C     | C        | C     | N          | N      |            |
|            | chi                | C     | C     | C        | C     | N          | N      |            |
|            | chi square         | C     | C     | C        | C     | N          | N      |            |
|            | cosine             | C     | C     | C        | T     | N          | N      |            |
|            | exponential        | C     | C     | C        | T     | N          | N      |            |
|            | dagum              | N     | N     | N        | N     | N          | N      | B          |
|            | davis              | N     | N     | N        | N     | N          | N      | B          |
|            | gamma              | C     | C     | C        | T     | N          | N      |            |
|            | gumbel             | C     | C     | C        | T     | N          | N      |            |
|            | gompertz           | C     | C     | C        | T     | N          | N      |            |
|            | fisher z           | C     | C     | C        | C     | N          | N      |            |
|            | fisk               | C     | C     | C        | C     | N          | N      |            |
|            | frechet            | C     | C     | C        | T     | N          | N      |            |
|            | half normal        | C     | C     | C        | T     | N          | N      |            |
|            | hyperbolic secant  | C     | C     | C        | T     | N          | N      |            |
|            | inverse gamma      | C     | C     | C        | T     | N          | N      |            |
|            | inverse gauss      | C     | C     | N        | T     | N          | N      |            |
|            | inverse chi        | C     | C     | C        | T     | N          | N      |            |
|            | inverse chi sq     | C     | C     | C        | T     | N          | N      |            |
|            | irwin hall         | C     | C     | C        | C     | N          | N      |            |
|            | johnson sb         | C     | C     | C        | T     | N          | N      |            |
|            | johnson su         | C     | C     | C        | T     | N          | N      |            |
|            | kumaraswamy        | C     | C     | C        | T     | N          | N      |            |
|            | laplace            | C     | C     | C        | T     | N          | N      |            |
|            | log logistic       | C     | C     | C        | T     | N          | N      |            |
|            | log normal         | C     | C     | C        | T     | N          | N      |            |
|            | logistic           | C     | C     | C        | T     | N          | N      |            |
|            | lomax              | C     | C     | C        | T     | N          | N      |            |
|            | maxwell            | C     | C     | C        | T     | N          | N      |            |
|            | nakagami           | C     | C     | C        | T     | N          | N      |            |
|            | non central chi sq | N     | N     | N        | N     | N          | N      | AAA        |
|            | non central f      | N     | N     | N        | N     | N          | N      | AAA        |
|            | non central t      | N     | N     | N        | N     | N          | N      | AAA        |
|            | pareto             | C     | C     | C        | T     | N          | N      |            |
|            | power              | C     | C     | C        | T     | N          | N      |            |
|            | rayleigh           | C     | C     | C        | T     | N          | N      |            |
|            | reciprocal         | C     | C     | C        | C     | N          | N      |            |
|            | rice               | N     | N     | N        | N     | N          | N      | AAA        |
|            | skew normal        | C     | C     | N        | T     | N          | N      |            |
|            | snedecor f         | C     | C     | C        | T     | N          | N      |            |
|            | student t          | C     | C     | C        | T     | N          | N      |            |
|            | suzuki             | N     | N     | N        | N     | N          | N      | B          |
|            | triangular         | C     | C     | C        | T     | N          | N      |            |
|            | uniform            | C     | C     | C        | T     | N          | N      |            |
|            | u quadratic        | C     | C     | C        | T     | N          | N      |            |
|            | voigt              | C     | C     | N        | T     | N          | N      |            |
|            | weibull            | C     | C     | C        | T     | N          | N      |            |
|            | wigner semicircle  | C     | C     | C        | T     | N          | N      |            |
|            | ***coverage***     | 48/56 | 48/56 | 45/56    | 47/56 | 0/56       | 0/56   |            |
| stable     | cauchy             | C     | C     | C        | T     | N          | N      |            |
|            | delta              | C     | C     | C        | T     | -          | N      |            |
|            | holtsmark          | C     | C     | C        | C     | N          | N      |            |
|            | landau             | C     | C     | C        | C     | N          | N      |            |
|            | map-airy           | C     | C     | C        | C     | N          | N      |            |
|            | levy               | C     | C     | C        | T     | N          | N      |            |
|            | normal             | C     | C     | C        | T     | N          | N      |            |
|            | sas point5         | C	  | C     | C        | C     | N          | N      |            |
|            | ***coverage***     | 8/8   | 8/8   | 8/8      | 8/8   | 0/8        | 0/8    |            |