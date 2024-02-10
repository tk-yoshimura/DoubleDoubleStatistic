using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ExponentialDistribution : ContinuousDistribution {

        public ddouble Lambda { get; }

        public ExponentialDistribution(ddouble lambda) {
            ValidateScale(lambda);

            Lambda = lambda;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble pdf = Lambda * Exp(-Lambda * x);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = -Expm1(-Lambda * x);

                return cdf;
            }
            else {
                ddouble cdf = Exp(-Lambda * x);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble quantile = -Log1p(-p) / Lambda;

                return quantile;
            }
            else {
                ddouble quantile = -Log(p) / Lambda;

                if (IsNegative(quantile)) {
                    return 0d;
                }

                return quantile;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => 1d / Lambda;
        public override ddouble Median => Ln2 / Lambda;
        public override ddouble Mode => 0d;

        public override ddouble Variance => 1d / (Lambda * Lambda);
        public override ddouble Skewness => 2;
        public override ddouble Kurtosis => 6;

        public override ddouble Entropy => 1 - Log(Lambda);

        public override string ToString() {
            return $"{typeof(ExponentialDistribution).Name}[lambda={Lambda}]";
        }
    }
}
