using DoubleDouble;
using DoubleDoubleRootFinding;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ChiSquaredDistribution : Distribution {

        public int K { get; }

        private readonly ddouble pdf_norm, cdf_norm;

        public ChiSquaredDistribution(int k) {
            ValidateShape(k, k => k > 0);

            this.K = k;

            this.cdf_norm = 1d / Gamma(k * 0.5d);
            this.pdf_norm = Pow2(k * -0.5d) * cdf_norm;
        }

        public override ddouble PDF(ddouble x) {
            if (IsZero(x)) {
                return K <= 1 ? PositiveInfinity : (K <= 2 ? Point5 : Zero);
            }

            ddouble pdf = pdf_norm * Exp((K * 0.5d - 1d) * Log(x) - Ldexp(x, -1));

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            ddouble cdf = cdf_norm * LowerIncompleteGamma(K * 0.5d, Ldexp(x, -1));

            cdf = IsFinite(cdf) ? cdf : (x > One ? One : Zero);

            return cdf;
        }

        public override (ddouble min, ddouble max) Support => (Zero, PositiveInfinity);

        public override ddouble Mean => K;
        public override ddouble Median {
            get {
                ddouble x0 = K * Cube(1d - (ddouble)2 / ((ddouble)9 * K));
                ddouble x = NewtonRaphsonFinder.RootFind(x => (CDF(x) - Point5, PDF(x)), x0, iters: 256);

                return x;
            }
        }
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
