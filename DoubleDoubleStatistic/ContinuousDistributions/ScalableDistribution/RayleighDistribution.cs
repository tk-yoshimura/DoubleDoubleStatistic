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
    public class RayleighDistribution : ScalableDistribution<RayleighDistribution>,
        IMultiplyOperators<RayleighDistribution, ddouble, RayleighDistribution>,
        IDivisionOperators<RayleighDistribution, ddouble, RayleighDistribution>,
        IFittableContinuousDistribution<RayleighDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv, sigma_sq;

        public RayleighDistribution() : this(sigma: 1d) { }

        public RayleighDistribution(ddouble sigma) {
            ValidateScale(sigma);

            Sigma = sigma;
            sigma_inv = 1d / sigma;
            sigma_sq = sigma * sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble pdf = u * Exp(u * u * -0.5d) * sigma_inv;
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
                if (u <= 0d || u2 <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(u2)) {
                    return 1d;
                }

                ddouble cdf = -Expm1(u2 * -0.5d);

                if (IsNaN(cdf)) {
                    return 1d;
                }

                return cdf;
            }
            else {
                if (u <= 0d || u2 <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(u2)) {
                    return 0d;
                }

                ddouble cdf = Exp(u2 * -0.5d);

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
                ddouble v = -2d * Log1p(-p);
                if (IsNegative(v)) {
                    return 0d;
                }

                ddouble x = Sigma * Sqrt(-2d * Log1p(-p));

                return x;
            }
            else {
                ddouble v = -2d * Log(p);
                if (IsNegative(v)) {
                    return 0d;
                }

                ddouble x = Sigma * Sqrt(v);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = (double)Sigma * double.Sqrt(-2d * double.Log(u));

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Sigma * Sqrt(PI * 0.5d);

        public override ddouble Median =>
            Sigma * Sqrt(2d * Ln2);

        public override ddouble Mode => Sigma;

        public override ddouble Variance =>
            (4d - PI) / 2d * sigma_sq;

        public override ddouble Skewness =>
            2d * Sqrt(PI) * (PI - 3d) / ExMath.Pow3d2(4d - PI);

        public override ddouble Kurtosis =>
            (-16d + PI * (24d + PI * -6d)) / Square(4d - PI);

        public override ddouble Entropy =>
            1d + Log(Sigma / Sqrt2) + EulerGamma * 0.5d;

        public static RayleighDistribution operator *(RayleighDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public static RayleighDistribution operator /(RayleighDistribution dist, ddouble k) {
            return new(dist.Sigma / k);
        }

        public static (RayleighDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (RayleighDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileScaleFitter<RayleighDistribution>.Fit(new RayleighDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(RayleighDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := u * Exp(-u^2 / 2) / sigma, u = x / sigma";
    }
}
