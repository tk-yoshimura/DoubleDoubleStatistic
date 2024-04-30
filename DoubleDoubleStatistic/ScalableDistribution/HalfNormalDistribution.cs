using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class HalfNormalDistribution : ScalableDistribution<HalfNormalDistribution>,
        IMultiplyOperators<HalfNormalDistribution, ddouble, HalfNormalDistribution>,
        IDivisionOperators<HalfNormalDistribution, ddouble, HalfNormalDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        public HalfNormalDistribution() : this(sigma: 1d) { }

        public HalfNormalDistribution(ddouble sigma) {
            ValidateScale(sigma);

            Sigma = sigma;

            sigma_sq = sigma * sigma;
            pdf_norm = Sqrt2 / (sigma * Sqrt(PI));
            exp_scale = -1d / (2d * sigma_sq);
            erf_scale = 1d / (Sqrt2 * sigma);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Exp(Square(x) * exp_scale);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

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

        public override ddouble Skewness => Sqrt2 * (4d - PI) / ExMath.Pow3d2(PI - 2d);

        public override ddouble Kurtosis => 8d * (PI - 3d) / Square(PI - 2d);

        public override ddouble Entropy => (Log(PI / 2) + 1d) * 0.5d + Log(Sigma);

        public static HalfNormalDistribution operator *(HalfNormalDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public static HalfNormalDistribution operator /(HalfNormalDistribution dist, ddouble k) {
            return new(dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(HalfNormalDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := exp(-u^2) * sqrt(2 / pi) / sigma, u = x / sigma, (x >= 0)";
    }
}
