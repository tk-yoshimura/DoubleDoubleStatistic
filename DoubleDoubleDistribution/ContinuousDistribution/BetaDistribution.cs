using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class BetaDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble pdf_norm;

        public BetaDistribution(ddouble alpha, ddouble beta) {
            ValidateShape(alpha, alpha => alpha > 0);
            ValidateShape(beta, beta => beta > 0);

            Alpha = alpha;
            Beta = beta;

            pdf_norm = 1d / Beta(alpha, beta);
        }

        public override ddouble PDF(ddouble x) {
            if (x < 0d || x > 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Pow(x, Alpha - 1d) * Beta(1d - x, Beta - 1d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = IncompleteBetaRegularized(x, Alpha, Beta);

                return cdf;
            }
            else {
                ddouble cdf = IncompleteBetaRegularized(1d - x, Beta, Alpha);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseIncompleteBeta(p, Alpha, Beta);

                return x;
            }
            else {
                ddouble x = InverseIncompleteBeta(1d - p, Alpha, Beta);

                return x;
            }
        }

        public override bool Symmetric => Alpha == Beta;

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean => Alpha / (Alpha + Beta);

        public override ddouble Median => InverseIncompleteBeta(0.5, Alpha, Beta);

        public override ddouble Mode =>
            Alpha > 1d && Beta > 1d ? (Alpha - 1d) / (Alpha + Beta - 2d) :
            Alpha <= 1d && Beta >= 1d && Alpha != Beta ? 0d :
            Alpha >= 1d && Beta <= 1d && Alpha != Beta ? 1d :
            NaN;

        public override ddouble Variance =>
            Alpha * Beta / (Square(Alpha + Beta) * (Alpha + Beta + 1d));

        public override ddouble Skewness =>
            2 * (Beta - Alpha) * Sqrt(Alpha + Beta + 1d) / ((Alpha + Beta + 2d) * Sqrt(Alpha * Beta));

        public override ddouble Kurtosis =>
            6 * (Square(Alpha - Beta) * (Alpha + Beta + 1d) - Alpha * Beta * (Alpha + Beta + 2d)) / (Alpha * Beta * (Alpha + Beta + 2d) * (Alpha + Beta + 3d));

        public override ddouble Entropy =>
            LogBeta(Alpha, Beta) - (Alpha - 1d) * Digamma(Alpha) - (Beta - 1d) * Digamma(Beta) + (Alpha + Beta - 2d) * Digamma(Alpha + Beta);

        public override string ToString() {
            return $"{typeof(BetaDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }
    }
}
