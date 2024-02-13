using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ExtremeValueDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        public ExtremeValueDistribution(ddouble alpha, ddouble beta) {
            ValidateLocation(alpha);
            ValidateShape(beta, beta => beta > 0);

            Alpha = alpha;
            Beta = beta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (Alpha - x) / Beta;

            ddouble pdf = Exp(-Exp(u) + u) / Beta;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (Alpha - x) / Beta;

            if (interval == Interval.Lower) {
                ddouble cdf = Exp(-Exp(u));

                return cdf;
            }
            else {
                ddouble cdf = -Expm1(-Exp(u));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = Log(-Log(p));
                ddouble quantile = Alpha - Beta * u;

                quantile = IsNaN(quantile) ? PositiveInfinity : quantile;

                return quantile;
            }
            else {
                ddouble u = Log(-Log1p(-p));
                ddouble quantile = Alpha - Beta * u;

                return quantile;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Alpha + EulerGamma * Beta;
        public override ddouble Median => Quantile(0.5);
        public override ddouble Mode => Alpha;

        public override ddouble Variance => Square(PI * Beta) / 6d;
        public override ddouble Skewness => 12d * Sqrt(6d) * Zeta3 / Cube(PI);
        public override ddouble Kurtosis => (ddouble)12 / 5;

        public override ddouble Entropy => Log(Beta) + EulerGamma + 1d;

        public override string ToString() {
            return $"{typeof(ExtremeValueDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }
    }
}
