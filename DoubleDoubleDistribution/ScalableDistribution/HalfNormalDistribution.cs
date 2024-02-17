using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class HalfNormalDistribution : ScalableDistribution<HalfNormalDistribution>,
        IMultiplyOperators<HalfNormalDistribution, ddouble, HalfNormalDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        public HalfNormalDistribution() : this(sigma: 1) { }

        public HalfNormalDistribution(ddouble sigma) {
            ValidateScale(sigma);

            Sigma = sigma;

            sigma_sq = sigma * sigma;
            pdf_norm = Sqrt2 / (sigma * Sqrt(PI));
            exp_scale = -1d / (2 * sigma_sq);
            erf_scale = 1d / (Sqrt2 * sigma);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Exp(Square(x) * exp_scale);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }

                ddouble cdf = Erf(x * erf_scale);

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                ddouble cdf = Erfc(x * erf_scale);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sigma * Sqrt2 * InverseErf(p);

                return x;
            }
            else {
                ddouble x = Sigma * Sqrt2 * InverseErfc(p);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Sigma * Sqrt2 / Sqrt(PI);

        public override ddouble Median => Sigma * Sqrt2 * InverseErf(0.5d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => sigma_sq * (1d - 2 * RcpPI);

        public override ddouble Skewness => Sqrt2 * (4d - PI) / Cube(Sqrt(PI - 2d));

        public override ddouble Kurtosis => 8 * (PI - 3d) / Square(PI - 2d);

        public override ddouble Entropy => Log2(2 * PI * E * sigma_sq) / 2 - 1d;

        public static HalfNormalDistribution operator *(HalfNormalDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public override string ToString() {
            return $"{typeof(HalfNormalDistribution).Name}[sigma={Sigma}]";
        }
    }
}
