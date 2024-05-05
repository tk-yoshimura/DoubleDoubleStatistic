using DoubleDouble;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BetaDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble pdf_norm;

        public BetaDistribution(ddouble alpha, ddouble beta) {
            ValidateShape(alpha, alpha => alpha > 0d);
            ValidateShape(beta, beta => beta > 0d);

            Alpha = alpha;
            Beta = beta;

            pdf_norm = 1d / Beta(alpha, beta);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x < 0d || x > 1d) {
                return 0d;
            }

            if (Alpha != 1d && Beta != 1d) {
                ddouble pdf = pdf_norm * Pow(x, Alpha - 1d) * Pow(1d - x, Beta - 1d);

                return pdf;
            }
            else if (Alpha != 1d) {
                ddouble pdf = pdf_norm * Pow(x, Alpha - 1d);

                return pdf;
            }
            else if (Beta != 1d) {
                ddouble pdf = pdf_norm * Pow(1d - x, Beta - 1d);

                return pdf;
            }
            else {
                return 1d;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (x >= 1d) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(x, Alpha, Beta);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (x >= 1d) {
                    return 0d;
                }

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

        public override ddouble Median => InverseIncompleteBeta(0.5d, Alpha, Beta);

        public override ddouble Mode =>
            Alpha > 1d && Beta > 1d ? (Alpha - 1d) / (Alpha + Beta - 2d) :
            Alpha <= 1d && Beta >= 1d && Alpha != Beta ? 0d :
            Alpha >= 1d && Beta <= 1d && Alpha != Beta ? 1d :
            NaN;

        public override ddouble Variance =>
            Alpha * Beta / (Square(Alpha + Beta) * (Alpha + Beta + 1d));

        public override ddouble Skewness =>
            2d * (Beta - Alpha) * Sqrt(Alpha + Beta + 1d) / ((Alpha + Beta + 2d) * Sqrt(Alpha * Beta));

        public override ddouble Kurtosis =>
            6d * (Square(Alpha - Beta) * (Alpha + Beta + 1d) - Alpha * Beta * (Alpha + Beta + 2d)) / (Alpha * Beta * (Alpha + Beta + 2d) * (Alpha + Beta + 3d));

        public override ddouble Entropy =>
            LogBeta(Alpha, Beta) - (Alpha - 1d) * Digamma(Alpha) - (Beta - 1d) * Digamma(Beta) + (Alpha + Beta - 2d) * Digamma(Alpha + Beta);

        public override string ToString() {
            return $"{typeof(BetaDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := x^(alpha - 1) * (1 - x)^(beta - 1) / beta(alpha, beta)";
    }
}
