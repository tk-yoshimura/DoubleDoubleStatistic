using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ChiSquareDistribution : ContinuousDistribution,
        IAdditionOperators<ChiSquareDistribution, ChiSquareDistribution, ChiSquareDistribution> {

        public int K { get; }

        private readonly ddouble c, pdf_lognorm;

        public ChiSquareDistribution(int k) {
            ValidateShape(k, k => k > 0);

            K = k;

            c = K * 0.5d - 1d;
            pdf_lognorm = k * 0.5d + LogGamma(k * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsZero(x)) {
                return K <= 1 ? PositiveInfinity : K <= 2 ? 0.5d : 0d;
            }

            ddouble pdf = Pow2(c * Log2(x) - Ldexp(x, -1) * LbE - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(K * 0.5d, Ldexp(x, -1));

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(K * 0.5d, Ldexp(x, -1));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseLowerIncompleteGamma(K * 0.5d, p) * 2;

                return x;
            }
            else {
                ddouble x = InverseUpperIncompleteGamma(K * 0.5d, p) * 2;

                return x;
            }
        }

        public override bool AdditiveClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => K;

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => Max(K - 2, 0);

        public override ddouble Variance => Ldexp(K, 1);

        public override ddouble Skewness => Sqrt((ddouble)8 / K);

        public override ddouble Kurtosis => (ddouble)12 / K;

        public override ddouble Entropy {
            get {
                ddouble k_half = Ldexp(K, -1);

                return k_half + Log(2 * Gamma(k_half)) + (1 - k_half) * Digamma(k_half);
            }
        }

        public static ChiSquareDistribution operator +(ChiSquareDistribution dist1, ChiSquareDistribution dist2) {
            return new(dist1.K + dist2.K);
        }

        public override string ToString() {
            return $"{typeof(ChiSquareDistribution).Name}[k={K}]";
        }
    }
}
