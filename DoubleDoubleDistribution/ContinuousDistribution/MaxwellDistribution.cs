using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class MaxwellDistribution : ContinuousDistribution,
        IMultiplyOperators<MaxwellDistribution, ddouble, MaxwellDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble sigma_sq, sigma_cb, pdf_norm;

        public MaxwellDistribution() : this(sigma: 1d) { }

        public MaxwellDistribution(ddouble sigma) {
            ValidateScale(sigma);

            Sigma = sigma;
            sigma_sq = sigma * sigma;
            sigma_cb = sigma * sigma * sigma;
            pdf_norm = Sqrt(2 * RcpPI) / sigma_cb;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble pdf = pdf_norm * x * x * Exp(-x * x / (2 * sigma_sq));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble x2 = x * x;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x2)) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(1.5d, Ldexp(x2 / sigma_sq, -1));

                if (IsNaN(cdf)) {
                    return 1d;
                }

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x2)) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(1.5d, Ldexp(x2 / sigma_sq, -1));

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
                ddouble x = Sqrt(InverseLowerIncompleteGamma(1.5d, p) * 2) * Sigma;

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(1.5d, p) * 2) * Sigma;

                return x;
            }
        }

        public override bool Scalable => true;

        public override (ddouble min, ddouble max) Support => (Zero, PositiveInfinity);

        public override ddouble Mean => 2 * Sigma * Sqrt(2 * RcpPI);
        public override ddouble Median => Quantile(0.5);
        public override ddouble Mode => Sqrt2 * Sigma;
        public override ddouble Variance => sigma_sq * (3d - 8d * RcpPI);
        public override ddouble Skewness => 2 * Sqrt2 * (16d - PI * 5) / Cube(Sqrt(3d * PI - 8d));

        public override ddouble Kurtosis => (-392d + PI * (160d + PI * -12d)) / Square(3d * PI - 8d);

        public override ddouble Entropy => -0.5d + Log(Sigma * Sqrt2 * Sqrt(PI)) + EulerGamma;

        public static MaxwellDistribution operator *(MaxwellDistribution dist, ddouble k) {
            return new(k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(MaxwellDistribution).Name}[sigma={Sigma}]";
        }
    }
}
