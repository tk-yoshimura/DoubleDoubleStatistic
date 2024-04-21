# DoubleDoubleDistribution
 Double-Double Continuous Distribution Implements

## Dev Stage
N: Not Implemented  
C: Test Checked  
T: TODO Test  

| category   | distribution       | pdf   | cdf   | quantile | stats | random gen | define | difficulty |
| ---------- | ------------------ | ----- | ----- | -------- | ----- | ---------- | ------ | ---------- |
| continuous | alpha              | C     | C     | C        | T     | N          | N      |            |
|            | arcsine            | C     | C     | C        | T     | N          | N      |            |
|            | argus              | N     | N     | N        | N     | N          | N      | A          |
|            | benktander         | N     | N     | N        | N     | N          | N      | C          |
|            | birnbaum saunders  | T     | T     | T        | T     | N          | N      |            |
|            | beta               | C     | C     | C        | T     | N          | N      |            |
|            | beta prime         | C     | C     | C        | T     | N          | N      |            |
|            | bradford           | C     | C     | C        | T     | N          | N      |            |
|            | burr               | C     | C     | C        | T     | N          | N      |            |
|            | chi                | C     | C     | C        | T     | N          | N      |            |
|            | chi square         | C     | C     | C        | T     | N          | N      |            |
|            | cosine             | C     | C     | C        | T     | N          | N      |            |
|            | exponential        | C     | C     | C        | T     | N          | N      |            |
|            | dagum              | N     | N     | N        | N     | N          | N      | B          |
|            | davis              | N     | N     | N        | N     | N          | N      | B          |
|            | gamma              | T     | T     | T        | T     | N          | N      |            |
|            | gumbel             | C     | C     | C        | T     | N          | N      |            |
|            | gompertz           | T     | T     | T        | T     | N          | N      |            |
|            | fisher z           | T     | T     | T        | T     | N          | N      |            |
|            | fisk               | C     | C     | C        | T     | N          | N      |            |
|            | frechet            | T     | T     | T        | T     | N          | N      |            |
|            | half normal        | T     | T     | T        | T     | N          | N      |            |
|            | hyperbolic secant  | C     | C     | C        | T     | N          | N      |            |
|            | inverse gamma      | T     | T     | N        | T     | N          | N      |            |
|            | inverse gauss      | T     | T     | N        | T     | N          | N      |            |
|            | inverse chi        | T     | T     | T        | T     | N          | N      |            |
|            | inverse chi sq     | T     | T     | T        | T     | N          | N      |            |
|            | irwin hall         | C     | C     | C        | C     | N          | N      |            |
|            | johnson sb         | T     | T     | T        | N     | N          | N      |            |
|            | johnson su         | T     | T     | T        | T     | N          | N      |            |
|            | kumaraswamy        | T     | T     | T        | T     | N          | N      |            |
|            | laplace            | C     | C     | C        | T     | N          | N      |            |
|            | log logistic       | T     | T     | T        | T     | N          | N      |            |
|            | log normal         | C     | C     | C        | T     | N          | N      |            |
|            | logistic           | C     | C     | C        | T     | N          | N      |            |
|            | lomax              | T     | T     | T        | T     | N          | N      |            |
|            | maxwell            | T     | T     | T        | T     | N          | N      |            |
|            | nakagami           | T     | T     | T        | T     | N          | N      |            |
|            | non central chi sq | N     | N     | N        | N     | N          | N      | AAA        |
|            | non central f      | N     | N     | N        | N     | N          | N      | AAA        |
|            | non central t      | N     | N     | N        | N     | N          | N      | AAA        |
|            | pareto             | T     | T     | T        | T     | N          | N      |            |
|            | power              | T     | T     | T        | T     | N          | N      |            |
|            | rayleigh           | T     | T     | T        | T     | N          | N      |            |
|            | reciprocal         | C     | C     | C        | C     | N          | N      |            |
|            | rice               | N     | N     | N        | N     | N          | N      | AAA        |
|            | skew normal        | C     | C     | N        | T     | N          | N      |            |
|            | snedecor f         | T     | T     | T        | T     | N          | N      |            |
|            | student t          | C     | C     | C        | T     | N          | N      |            |
|            | suzuki             | N     | N     | N        | N     | N          | N      | B          |
|            | triangular         | T     | T     | T        | T     | N          | N      |            |
|            | uniform            | T     | T     | T        | T     | N          | N      |            |
|            | u shape            | T     | T     | T        | T     | N          | N      |            |
|            | voigt              | C     | C     | N        | T     | N          | N      |            |
|            | weibull            | T     | T     | T        | T     | N          | N      |            |
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