using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class GammaDistribution : ContinuousDistribution {

        public ddouble K { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_lognorm;

        public GammaDistribution(ddouble k, ddouble theta) {
            ValidateShape(k, k => k > 0);
            ValidateScale(theta);

            K = k;
            Theta = theta;

            pdf_lognorm = -LogGamma(K) + K * Log(Theta);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble pdf = Exp((K - 1d) * Log(x) - x / Theta + pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = LowerIncompleteGammaRegularized(K, x / Theta);

                return cdf;
            }
            else {
                ddouble cdf = UpperIncompleteGammaRegularized(K, x / Theta);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble quantile = InverseLowerIncompleteGamma(K, p) * Theta;

                return quantile;
            }
            else {
                ddouble quantile = InverseUpperIncompleteGamma(K, p) * Theta;

                return quantile;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => K * Theta;
        public override ddouble Median => InverseLowerIncompleteGamma(K, 0.5) * Theta;
        public override ddouble Mode => K >= 1d ? (K - 1d) * Theta : 0d;

        public override ddouble Variance => K * Theta * Theta;
        public override ddouble Skewness => 2 / Sqrt(K);
        public override ddouble Kurtosis => 6 / K;

        public override ddouble Entropy =>
            K + Log(Theta) + LogGamma(K) - (K - 1d) * Digamma(K);

        public override string ToString() {
            return $"{typeof(GammaDistribution).Name}[k={K},theta={Theta}]";
        }
    }
}
