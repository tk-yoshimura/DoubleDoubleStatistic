using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class GammaDistribution : ScalableDistribution<GammaDistribution>,
        IMultiplyOperators<GammaDistribution, ddouble, GammaDistribution> {

        public ddouble K { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_lognorm, theta_inv;

        public GammaDistribution(ddouble k, ddouble theta) {
            ValidateShape(k, k => k > 0);
            ValidateScale(theta);

            K = k;
            Theta = theta;

            pdf_lognorm = -LogGamma(K) + K * Log(Theta);
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble pdf = Exp((K - 1d) * Log(x) - x * theta_inv + pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = LowerIncompleteGammaRegularized(K, x * theta_inv);

                return cdf;
            }
            else {
                ddouble cdf = UpperIncompleteGammaRegularized(K, x * theta_inv);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseLowerIncompleteGamma(K, p) * Theta;

                return x;
            }
            else {
                ddouble x = InverseUpperIncompleteGamma(K, p) * Theta;

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => K * Theta;

        public override ddouble Median =>
            InverseLowerIncompleteGamma(K, 0.5) * Theta;

        public override ddouble Mode => K >= 1d ? (K - 1d) * Theta : 0d;

        public override ddouble Variance => K * Theta * Theta;

        public override ddouble Skewness => 2 / Sqrt(K);

        public override ddouble Kurtosis => 6 / K;

        public override ddouble Entropy =>
            K + Log(Theta) + LogGamma(K) - (K - 1d) * Digamma(K);

        public static GammaDistribution operator *(GammaDistribution dist, ddouble k) {
            return new(dist.K, dist.Theta * k);
        }

        public override string ToString() {
            return $"{typeof(GammaDistribution).Name}[k={K},theta={Theta}]";
        }
    }
}
