using DoubleDouble;
using DoubleDoubleRootFinding;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ChiSquaredDistribution : Distribution {

        public int K { get; }

        private readonly ddouble c, pdf_norm;

        public ChiSquaredDistribution(int k) {
            ValidateShape(k, k => k > 0);

            this.K = k;

            this.c = K * 0.5d - 1d;
            this.pdf_norm = k * 0.5d + LogGamma(k * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsZero(x)) {
                return K <= 1 ? PositiveInfinity : (K <= 2 ? Point5 : Zero);
            }

            ddouble pdf = Pow2(c * Log2(x) - Ldexp(x, -1) * LbE - pdf_norm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = LowerIncompleteGammaRegularized(K * 0.5d, Ldexp(x, -1));

                return cdf;
            }
            else {
                ddouble cdf = UpperIncompleteGammaRegularized(K * 0.5d, Ldexp(x, -1));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble x = InverseLowerIncompleteGamma(K * 0.5d, p);

                return x;
            }
            else {
                ddouble x = InverseUpperIncompleteGamma(K * 0.5d, p);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (Zero, PositiveInfinity);

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

        public override string ToString() {
            return $"{typeof(ChiSquaredDistribution).Name}[k={K}]";
        }
    }
}
