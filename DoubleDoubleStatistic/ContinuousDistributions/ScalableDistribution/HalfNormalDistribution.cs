using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class HalfNormalDistribution : ScalableDistribution<HalfNormalDistribution>,
        IMultiplyOperators<HalfNormalDistribution, ddouble, HalfNormalDistribution>,
        IDivisionOperators<HalfNormalDistribution, ddouble, HalfNormalDistribution>,
        IFittableContinuousDistribution<HalfNormalDistribution> {

        public ddouble Sigma { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private readonly ddouble pdf_norm, sigma_inv;

        public HalfNormalDistribution() : this(sigma: 1d) { }

        public HalfNormalDistribution(ddouble sigma) {
            ValidateScale(sigma);

            Sigma = sigma;

            sigma_inv = 1d / sigma;
            pdf_norm = Sqrt2 / (sigma * Sqrt(PI));
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Exp(u * u * -0.5d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }

                ddouble cdf = Erf(u * sqrt2_inv);

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }

                ddouble cdf = Erfc(u * sqrt2_inv);

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

        public override double Sample(Random random) {
            double u = random.NextGaussian();

            double v = double.Abs(u) * (double)Sigma;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Sigma * Sqrt2 / Sqrt(PI);

        public override ddouble Median => Sigma * Sqrt2 * InverseErf(0.5d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => Sigma * Sigma * (1d - 2 * RcpPI);

        public override ddouble Skewness => Sqrt2 * (4d - PI) / ExMath.Pow3d2(PI - 2d);

        public override ddouble Kurtosis => 8d * (PI - 3d) / Square(PI - 2d);

        public override ddouble Entropy => (Log(PI / 2) + 1d) * 0.5d + Log(Sigma);

        public static HalfNormalDistribution operator *(HalfNormalDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public static HalfNormalDistribution operator /(HalfNormalDistribution dist, ddouble k) {
            return new(dist.Sigma / k);
        }

        public static (HalfNormalDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (HalfNormalDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileScaleFitter<HalfNormalDistribution>.Fit(new HalfNormalDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(HalfNormalDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := exp(-u^2) * sqrt(2 / pi) / sigma, u = x / sigma, (x >= 0)";
    }
}