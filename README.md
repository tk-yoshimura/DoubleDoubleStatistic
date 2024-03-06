# DoubleDoubleDistribution
 Double-Double Continuous Distribution Implements

## Dev Stage
N: Not Implemented  
C: Test Checked  
T: TODO Test  

| category   | distribution       | pdf   | cdf   | quantile | stats | random gen | define | difficulty |
| ---------- | ------------------ | ----- | ----- | -------- | ----- | ---------- | ------ | ---------- |
| continuous | alpha              | T     | T     | T        | T     | N          | N      |            |
|            | arcsine            | T     | T     | T        | T     | N          | N      |            |
|            | argus              | N     | N     | N        | N     | N          | N      | A          |
|            | birnbaum saunders  | N     | N     | N        | N     | N          | N      | C          |
|            | bradford           | T     | T     | T        | T     | N          | N      |            |
|            | beta               | T     | T     | T        | T     | N          | N      |            |
|            | beta prime         | T     | T     | T        | T     | N          | N      |            |
|            | burr               | T     | T     | T        | T     | N          | N      |            |
|            | chi                | T     | T     | T        | T     | N          | N      |            |
|            | chi square         | C     | C     | C        | T     | N          | N      |            |
|            | cosine             | T     | T     | T        | T     | N          | N      |            |
|            | exponential        | T     | T     | T        | T     | N          | N      |            |
|            | dagum              | N     | N     | N        | N     | N          | N      | B          |
|            | gamma              | T     | T     | T        | T     | N          | N      |            |
|            | gumbel             | T     | T     | T        | T     | N          | N      |            |
|            | gompertz           | T     | T     | T        | T     | N          | N      |            |
|            | fisher z           | T     | T     | T        | T     | N          | N      |            |
|            | fisk               | T     | T     | T        | T     | N          | N      |            |
|            | frechet            | T     | T     | T        | T     | N          | N      |            |
|            | half normal        | T     | T     | T        | T     | N          | N      |            |
|            | hyperbolic secant  | T     | T     | T        | T     | N          | N      |            |
|            | inverse gamma      | T     | T     | N        | T     | N          | N      |            |
|            | inverse gauss      | T     | T     | N        | T     | N          | N      |            |
|            | inverse chi        | T     | T     | T        | T     | N          | N      |            |
|            | inverse chi sq     | T     | T     | T        | T     | N          | N      |            |
|            | irwin hall         | T     | T     | T        | T     | N          | N      |            |
|            | johnson sb         | T     | T     | T        | N     | N          | N      |            |
|            | johnson su         | T     | T     | T        | T     | N          | N      |            |
|            | kumaraswamy        | T     | T     | T        | T     | N          | N      |            |
|            | laplace            | T     | T     | T        | T     | N          | N      |            |
|            | log logistic       | T     | T     | T        | T     | N          | N      |            |
|            | log normal         | C     | C     | C        | T     | N          | N      |            |
|            | logistic           | T     | T     | T        | T     | N          | N      |            |
|            | lomax              | T     | T     | T        | T     | N          | N      |            |
|            | maxwell            | T     | T     | T        | T     | N          | N      |            |
|            | nakagami           | T     | T     | T        | T     | N          | N      |            |
|            | non central chi sq | N     | N     | N        | N     | N          | N      | AAA        |
|            | non central f      | N     | N     | N        | N     | N          | N      | AAA        |
|            | non central t      | N     | N     | N        | N     | N          | N      | AAA        |
|            | pareto             | T     | T     | T        | T     | N          | N      |            |
|            | power              | T     | T     | T        | T     | N          | N      |            |
|            | rayleigh           | T     | T     | T        | T     | N          | N      |            |
|            | reciprocal         | T     | T     | T        | C     | N          | N      |            |
|            | rice               | N     | N     | N        | N     | N          | N      | AAA        |
|            | skew normal        | T     | T     | N        | T     | N          | N      |            |
|            | snedecor f         | T     | T     | T        | T     | N          | N      |            |
|            | student t          | C     | C     | C        | T     | N          | N      |            |
|            | triangular         | T     | T     | T        | T     | N          | N      |            |
|            | uniform            | T     | T     | T        | T     | N          | N      |            |
|            | u shape            | T     | T     | T        | T     | N          | N      |            |
|            | voigt              | C     | C     | N        | T     | N          | N      |            |
|            | weibull            | T     | T     | T        | T     | N          | N      |            |
|            | wigner semicircle  | C     | C     | C        | T     | N          | N      |            |
|            | ***coverage***     | 47/53 | 47/53 | 44/53    | 46/53 | 0/53       | 0/53   |            |
| stable     | cauchy             | C     | C     | C        | T     | N          | N      |            |
|            | delta              | C     | C     | C        | T     | -          | N      |            |
|            | holtsmark          | N     | N     | N        | N     | N          | N      | A          |
|            | landau             | N     | N     | N        | N     | N          | N      | A          |
|            | map-airy           | N     | N     | N        | N     | N          | N      | A          |
|            | levy               | T     | T     | T        | T     | N          | N      |            |
|            | normal             | C     | C     | C        | T     | N          | N      |            |
|            | sas point5         | N     | N     | N        | N     | N          | N      | A          |
|            | ***coverage***     | 4/8   | 4/8   | 4/8      | 4/8   | 0/8        | 0/8    |            |