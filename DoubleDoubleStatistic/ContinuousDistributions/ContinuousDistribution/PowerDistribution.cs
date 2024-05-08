using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class PowerDistribution : ContinuousDistribution {

        public ddouble K { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_lognorm, pdf_norm, k_inv, alpha_inv;

        public PowerDistribution(ddouble k, ddouble alpha) {
            ValidateShape(k, k => k > 0d);
            ValidateShape(alpha, alpha => alpha > 0d);

            K = k;
            Alpha = alpha;
            pdf_norm = alpha * Pow(K, alpha);
            pdf_lognorm = Log2(alpha) + Log2(k) * alpha;
            k_inv = 1d / k;
            alpha_inv = 1d / alpha;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || x > k_inv) {
                return 0d;
            }

            if (Alpha <= 2d) {
                ddouble pdf = pdf_norm * Pow(x, Alpha - 1d);

                return pdf;
            }
            else {
                ddouble pdf = Pow2(pdf_lognorm + Log2(x) * (Alpha - 1d));

                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                if (x > k_inv) {
                    return 1d;
                }

                ddouble cdf = Pow(K * x, Alpha);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                if (x > k_inv) {
                    return 0d;
                }

                ddouble cdf = 1d - Pow(K * x, Alpha);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p);
            }

            ddouble x = Pow(p, alpha_inv) * k_inv;

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(u, (double)alpha_inv) * (double)k_inv;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, k_inv);

        public override ddouble Mean => Alpha / (K * (Alpha + 1d));

        public override ddouble Median => 1d / (Pow2(alpha_inv) * K);

        public override ddouble Mode => Alpha < 1d ? 0d : Alpha > 1d ? k_inv : NaN;

        public override ddouble Variance =>
            Alpha / (Square(K * (Alpha + 1d)) * (Alpha + 2d));

        public override ddouble Skewness =>
            2d * (1d - Alpha) * Sqrt((Alpha + 2d) * alpha_inv) / (Alpha + 3d);

        public override ddouble Kurtosis =>
            6d * (2d + Alpha * (-6d + Alpha * (-1d + Alpha))) / (Alpha * (Alpha + 3d) * (Alpha + 4d));

        public override ddouble Entropy =>
            1d - Log(Alpha * K) - alpha_inv;

        public override string ToString() {
            return $"{typeof(PowerDistribution).Name}[k={K},alpha={Alpha}]";
        }

        public override string Formula => "p(x; k, alpha) := x^(alpha - 1) * alpha * k^alpha";
    }
}
