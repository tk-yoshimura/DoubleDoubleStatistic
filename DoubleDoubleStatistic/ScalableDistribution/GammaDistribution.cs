using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class GammaDistribution : ScalableDistribution<GammaDistribution>,
        IMultiplyOperators<GammaDistribution, ddouble, GammaDistribution>,
        IDivisionOperators<GammaDistribution, ddouble, GammaDistribution> {

        public ddouble K { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_lognorm, theta_inv;

        public GammaDistribution(ddouble k) : this(k, theta: 1d) { }

        public GammaDistribution(ddouble k, ddouble theta) {
            ValidateShape(k, k => k > 0d);
            ValidateScale(theta);

            K = k;
            Theta = theta;

            pdf_lognorm = LogGamma(k) * LbE;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u) || IsPositiveInfinity(u)) {
                return 0d;
            }

            if (IsZero(u)) {
                return K < 1d ? PositiveInfinity : K == 1d ? theta_inv : 0d;
            }

            ddouble pdf = Pow2((K - 1d) * Log2(u) - u * LbE - pdf_lognorm) * theta_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(K, u);

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(K, u);

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

        public static GammaDistribution operator /(GammaDistribution dist, ddouble k) {
            return new(dist.K, dist.Theta / k);
        }

        public override string ToString() {
            return $"{typeof(GammaDistribution).Name}[k={K},theta={Theta}]";
        }

        public override string Formula => "p(x; k, theta) := u^(k - 1) * exp(-u) / (gamma(k) * theta), u = x / theta";
    }
}
