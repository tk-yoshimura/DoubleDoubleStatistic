using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class PowerDistribution : ContinuousDistribution {

        public ddouble K { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_lognorm, k_inv, alpha_inv;

        public PowerDistribution(ddouble k, ddouble alpha) {
            ValidateShape(k, k => k > 0d);
            ValidateShape(alpha, alpha => alpha > 0d);

            K = k;
            Alpha = alpha;
            pdf_lognorm = Log(alpha) + Log(k) * alpha;
            k_inv = 1d / k;
            alpha_inv = 1d / alpha;
        }

        public override ddouble PDF(ddouble x) {
            if (x <= 0d || x > k_inv) {
                return 0d;
            }

            ddouble pdf = Exp(pdf_lognorm + Log(x) * (Alpha - 1d));
            
            return pdf;
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

        public override (ddouble min, ddouble max) Support => (0d, k_inv);

        public override ddouble Mean => Alpha / (K * (Alpha + 1d));

        public override ddouble Median => 1d / (Pow2(alpha_inv) * K);

        public override ddouble Mode => k_inv;

        public override ddouble Variance =>
            Alpha / (Square(K * (Alpha + 1d)) * (Alpha + 2d));

        public override ddouble Skewness =>
            2 * (1d - Alpha) * Sqrt((Alpha + 2) * alpha_inv) / (Alpha + 3d);

        public override ddouble Kurtosis =>
            6d * (2d + Alpha * (-6d + Alpha * (-1d + Alpha))) / (Alpha * (Alpha + 3d) * (Alpha + 4d));

        public override ddouble Entropy =>
            1d - Log(Alpha * K) - alpha_inv;

        public override string ToString() {
            return $"{typeof(PowerDistribution).Name}[k={K},alpha={Alpha}]";
        }
    }
}
