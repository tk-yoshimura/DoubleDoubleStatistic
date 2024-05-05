using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class KumaraswamyDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble ab, alpha_inv, beta_inv;

        public KumaraswamyDistribution(ddouble alpha, ddouble beta) {
            ValidateShape(alpha, alpha => alpha > 0d);
            ValidateShape(beta, beta => beta > 0d);

            Alpha = alpha;
            Beta = beta;

            ab = alpha * beta;
            alpha_inv = 1d / alpha;
            beta_inv = 1d / beta;
        }

        public override ddouble PDF(ddouble x) {
            if (x < 0d || x > 1d) {
                return 0d;
            }

            if (Beta == 1d) {
                return Alpha * Pow(x, Alpha - 1d);
            }
            if (Alpha == 1d) {
                ddouble xm1 = 1d - x;

                return Beta * Pow(xm1, Beta - 1d);
            }

            if (x == 0d) {
                return Alpha < 1d ? PositiveInfinity : 0d;
            }

            if (x == 1d) {
                return Beta < 1d ? PositiveInfinity : 0d;
            }

            ddouble xa = Pow(x, Alpha);

            ddouble pdf = ab * xa * Pow(1d - xa, Beta - 1d) / x;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }
                if (x > 1d) {
                    return 1d;
                }

                ddouble cdf = 1d - Pow(1d - Pow(x, Alpha), Beta);

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }
                if (x > 1d) {
                    return 0d;
                }

                ddouble cdf = Pow(1d - Pow(x, Alpha), Beta);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return Quantile(1d - p, Interval.Upper);
            }
            else {
                ddouble u = 1d - Pow(p, beta_inv);

                if (u <= 0d) {
                    return 0d;
                }
                if (u >= 1d) {
                    return 1d;
                }

                ddouble x = Pow(u, alpha_inv);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean => Beta * Beta(1d + alpha_inv, Beta);

        public override ddouble Median => Pow(1d - Pow2(-beta_inv), alpha_inv);

        public override ddouble Mode =>
            Alpha >= 1d && Beta >= 1d && !(Alpha == 1d && Beta == 1d)
            ? Pow((Alpha - 1d) / (ab - 1d), alpha_inv)
            : NaN;

        public override ddouble Variance =>
            Beta * Beta(Beta, 1d + 2d * alpha_inv) - Square(Mean);

        public override ddouble Skewness {
            get {
                ddouble b1 = Beta(Beta, 1d + alpha_inv);
                ddouble b2 = Beta(Beta, 1d + 2d * alpha_inv);
                ddouble b3 = Beta(Beta, 1d + 3d * alpha_inv);

                return (2d * Cube(Beta * b1) - 3d * Beta * Beta * b1 * b2 + Beta * b3)
                     / ExMath.Pow3d2(Beta * b2 - Square(Beta * b1));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble b1 = Beta(Beta, 1d + alpha_inv);
                ddouble b2 = Beta(Beta, 1d + 2d * alpha_inv);
                ddouble b3 = Beta(Beta, 1d + 3d * alpha_inv);
                ddouble b4 = Beta(Beta, 1d + 4d * alpha_inv);

                return (-3d * Pow(Beta * b1, 4) + 6d * Cube(Beta) * b1 * b1 * b2 - 4d * Beta * Beta * b1 * b3 + Beta * b4)
                    / Square(Beta * b2 - Square(Beta * b1)) - 3d;
            }
        }

        public override ddouble Entropy =>
            (1d - beta_inv) + (1d - alpha_inv) * (Digamma(Beta + 1d) + EulerGamma) - Log(ab);

        public override string ToString() {
            return $"{typeof(KumaraswamyDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := alpha * beta * x^(alpha - 1) * (1 - x^alpha)^(beta - 1)";
    }
}
