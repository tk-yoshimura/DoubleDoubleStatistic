using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class InverseGammaDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble pdf_lognorm;

        public InverseGammaDistribution() : this(alpha: 1d, beta: 1d) { }

        public InverseGammaDistribution(ddouble alpha, ddouble beta) {
            ValidateShape(alpha, alpha => alpha > 0d);
            ValidateShape(beta, beta => beta > 0d);

            Alpha = alpha;
            Beta = beta;

            pdf_lognorm = Alpha * Log2(Beta) - LogGamma(Alpha) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (x <= 0d) {
                return 0d;
            }

            ddouble pdf = Pow2(-Beta / x * LbE - (Alpha + 1d) * Log2(x) + pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Alpha, Beta / x);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Alpha, Beta / x);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return 0d;
                }

                ddouble x = Beta / InverseUpperIncompleteGamma(Alpha, p);

                return x;
            }
            else {
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Beta / InverseLowerIncompleteGamma(Alpha, p);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d)
            ? (Beta / (Alpha - 1d))
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode =>
            Beta / (Alpha + 1d);

        public override ddouble Variance => (Alpha > 2d)
            ? Square(Beta / (Alpha - 1d)) / (Alpha - 2d)
            : NaN;

        public override ddouble Skewness => (Alpha > 3d)
            ? (4 * Sqrt(Alpha - 2d) / (Alpha - 3d))
            : NaN;

        public override ddouble Kurtosis => (Alpha > 4d)
            ? ((30d * Alpha - 66d) / ((Alpha - 3d) * (Alpha - 4d)))
            : NaN;

        public override ddouble Entropy => Alpha + Log(Beta) + LogGamma(Alpha) - (Alpha + 1d) * Digamma(Alpha);

        public override string ToString() {
            return $"{typeof(InverseGammaDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := x^(- alpha - 1) * exp(-beta / x) * beta^alpha / gamma(alpha)";
    }
}
