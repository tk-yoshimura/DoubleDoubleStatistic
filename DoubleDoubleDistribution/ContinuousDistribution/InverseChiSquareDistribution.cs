using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class InverseChiSquareDistribution : ContinuousDistribution {

        public int K { get; }

        private readonly ddouble c, pdf_lognorm;

        public InverseChiSquareDistribution(int k) {
            ValidateShape(k, k => k > 0);

            K = k;

            c = K * 0.5d + 1d;
            pdf_lognorm = k * 0.5d + LogGamma(k * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble pdf = Pow2(-c * Log2(x) - LbE / (2 * x) - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(K * 0.5d, 1d / (2 * x));
            
                if (IsNaN(cdf)) {
                    return x < Mean ? 0d : 1d;
                }

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(K * 0.5d, 1d / (2 * x));

                if (IsNaN(cdf)) {
                    return x < Mean ? 1d : 0d;
                }

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = 0.5d / InverseUpperIncompleteGamma(K * 0.5d, p);

                return x;
            }
            else {
                ddouble x = 0.5d / InverseLowerIncompleteGamma(K * 0.5d, p);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => K > 2
            ? 1d / (ddouble)(K - 1)
            : NaN;

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => 1d / (ddouble)checked(K + 2);

        public override ddouble Variance => K > 4
            ? 2 / (Square(K - 2) * (ddouble)(K - 4))
            : NaN;

        public override ddouble Skewness => K > 6
            ? 4 * Sqrt(2 * (K - 4)) / (K - 6)
            : NaN;

        public override ddouble Kurtosis => K > 8
            ? (ddouble)checked(12 * (5 * K - 22)) / checked((K - 6) * (K - 8))
            : NaN;

        public override ddouble Entropy {
            get {
                ddouble k_half = Ldexp(K, -1);

                return k_half + Log(k_half * Gamma(k_half)) + (1 + k_half) * Digamma(k_half);
            }
        }

        public override string ToString() {
            return $"{typeof(InverseChiSquareDistribution).Name}[k={K}]";
        }
    }
}
