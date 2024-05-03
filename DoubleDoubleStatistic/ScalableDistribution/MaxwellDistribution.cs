using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class MaxwellDistribution : ScalableDistribution<MaxwellDistribution>,
        IMultiplyOperators<MaxwellDistribution, ddouble, MaxwellDistribution>,
        IDivisionOperators<MaxwellDistribution, ddouble, MaxwellDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv, sigma_sq, pdf_norm;

        public MaxwellDistribution() : this(sigma: 1d) { }

        public MaxwellDistribution(ddouble sigma) {
            ValidateScale(sigma);

            Sigma = sigma;
            sigma_inv = 1d / sigma;
            sigma_sq = sigma * sigma;
            pdf_norm = Sqrt(2d * RcpPI) * sigma_inv;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble u2 = u * u;

            ddouble pdf = pdf_norm * u2 * Exp(-u2 * 0.5d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;
            ddouble u2 = u * u;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(u2)) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(1.5d, u2 * 0.5d);

                if (IsNaN(cdf)) {
                    return 1d;
                }

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(u2)) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(1.5d, u2 * 0.5d);

                if (IsNaN(cdf)) {
                    return 0d;
                }

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(InverseLowerIncompleteGamma(1.5d, p) * 2d) * Sigma;

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(1.5d, p) * 2d) * Sigma;

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            2d * Sigma * Sqrt(2d * RcpPI);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Sqrt2 * Sigma;

        public override ddouble Variance =>
            sigma_sq * (3d - 8d * RcpPI);

        public override ddouble Skewness =>
            2 * Sqrt2 * (16d - PI * 5d) / ExMath.Pow3d2(3d * PI - 8d);

        public override ddouble Kurtosis =>
            (-384d + PI * (160d + PI * -12d)) / Square(3d * PI - 8d);

        public override ddouble Entropy =>
            -0.5d + Log(Sigma * Sqrt2 * Sqrt(PI)) + EulerGamma;

        public static MaxwellDistribution operator *(MaxwellDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public static MaxwellDistribution operator /(MaxwellDistribution dist, ddouble k) {
            return new(dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(MaxwellDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := sqrt(2 / pi) * u^2 * Exp(-u^2 / 2) / sigma, u = x / sigma";
    }
}
