using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class SnedecorFDistribution : ContinuousDistribution {

        public int D1 { get; }
        public int D2 { get; }

        private readonly ddouble pdf_norm, logd2xd2;

        public SnedecorFDistribution(int d1, int d2) {
            ValidateShape(d1, d1 => d1 > 0);
            ValidateShape(d2, d2 => d2 > 0);
            _ = checked(d1 + d2);

            D1 = d1;
            D2 = d2;

            pdf_norm = 1d / Beta(d1 * 0.5d, d2 * 0.5d);
            logd2xd2 = D2 * Log(D2);
        }

        public override ddouble PDF(ddouble x) {
            if (x <= 0d) {
                return 0d;
            }

            ddouble d1x = D1 * x;

            ddouble u = Exp(D1 * Log(d1x) + logd2xd2 - (D1 + D2) * Log(d1x + D2));
            ddouble pdf = this.pdf_norm * Sqrt(u) / x;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble d1x = D1 * x;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x)) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(d1x / (d1x + D2), D1 * 0.5d, D2 * 0.5d);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble cdf = IncompleteBetaRegularized(1d - d1x / (d1x + D2), D2 * 0.5d, D1 * 0.5d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = InverseIncompleteBeta(p, D1 * 0.5d, D2 * 0.5d);
                ddouble quantile = (D2 * u) / (D1 * (1d - u));

                return quantile;
            }
            else {
                ddouble u = InverseIncompleteBeta(1d - p, D1 * 0.5d, D2 * 0.5d);
                ddouble quantile = (D2 * u) / (D1 * (1d - u));

                return quantile;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => (D2 > 2) ? (ddouble)D2 / (D2 - 2) : NaN;
        public override ddouble Median => Quantile(0.5d);
        public override ddouble Mode => (D1 > 2) ? (ddouble)checked((D1 - 2) * D2) / (D1 * (D2 + 2)) : NaN;

        public override ddouble Variance => (D2 > 4)
            ? (ddouble)checked(2 * D2 * D2 * (D1 + D2 - 2)) / checked(D1 * (D2 - 2) * (D2 - 2) * (D2 - 4))
            : NaN;
        public override ddouble Skewness => (D2 > 6)
            ? checked(2 * D1 + D2 - 2) * Sqrt(8 * (D2 - 4)) / ((D2 - 6) * Sqrt(checked(D1 * (D1 + D2 - 2))))
            : NaN;
        public override ddouble Kurtosis => (D2 > 8)
            ? checked(12d * (ddouble)checked(D1 * (5 * D2 - 22) * (D1 + D2 - 2) + (D2 - 4) * (D2 - 2) * (D2 - 2))
                / checked(D1 * (D2 - 6) * (D2 - 8) * (D1 + D2 - 2)))
            : NaN;

        public override ddouble Entropy =>
            LogGamma(D1 * 0.5d) + LogGamma(D2 * 0.5d) - LogGamma((D1 + D2) * 0.5d)
            + (1d - D1 * 0.5d) * Digamma(1d + D1 * 0.5d) - (1d - D2 * 0.5d) * Digamma(1d + D2 * 0.5d)
            + ((D1 + D2) * 0.5d) * Digamma((D1 + D2) * 0.5d) + Log((ddouble)D1 / D2);

        public override string ToString() {
            return $"{typeof(SnedecorFDistribution).Name}[d1={D1},d2={D2}]";
        }
    }
}
