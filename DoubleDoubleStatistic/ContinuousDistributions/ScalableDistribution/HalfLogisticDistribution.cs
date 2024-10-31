using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class HalfLogisticDistribution : ScalableDistribution<HalfLogisticDistribution>,
        IMultiplyOperators<HalfLogisticDistribution, ddouble, HalfLogisticDistribution>,
        IDivisionOperators<HalfLogisticDistribution, ddouble, HalfLogisticDistribution>,
        IFittableContinuousDistribution<HalfLogisticDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public HalfLogisticDistribution() : this(sigma: 1d) { }

        public HalfLogisticDistribution(ddouble sigma) {
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Sigma = sigma;

            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv, v = Exp(-u);

            if (IsNaN(u)) {
                return NaN;
            }

            if (IsNegative(u) || IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble pdf = IsFinite(v) ? (2d * v) / Square(1d + v) * sigma_inv : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;
            ddouble v = Exp(-u);

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }
                if (IsPositiveInfinity(u)) {
                    return 1d;
                }


                ddouble cdf = (1d - v) / (1d + v);

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = (2d * v) / (1d + v);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p >= 1d) {
                    return PositiveInfinity;
                }

                ddouble u = Log((1d + p) / (1d - p));

                ddouble x = u * Sigma;

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                ddouble u = Log((2d - p) / p);

                ddouble x = u * Sigma;

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Log((2d - u) / u) * (double)Sigma;

            return v;
        }

        public override ddouble Mean => Log(4d) * Sigma;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => (Pi * Pi / 3d - Square(Log(4d))) * Sigma * Sigma;

        public override ddouble Skewness => "1.540328803404880246003799149979";

        public override ddouble Kurtosis => "3.583735664456714756880688027286";

        public override ddouble Entropy => Log(Sigma / 2d) + 2d;

        public static HalfLogisticDistribution operator *(HalfLogisticDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public static HalfLogisticDistribution operator /(HalfLogisticDistribution dist, ddouble k) {
            return new(dist.Sigma / k);
        }

        public static (HalfLogisticDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (HalfLogisticDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileScaleFitter<HalfLogisticDistribution>.Fit(new HalfLogisticDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(HalfLogisticDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := 2 * exp(-u) / (1 + exp(-u))^2 / sigma, u = x / sigma";
    }
}
