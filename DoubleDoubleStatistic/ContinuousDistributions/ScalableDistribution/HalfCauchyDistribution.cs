using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class HalfCauchyDistribution : ScalableDistribution<HalfCauchyDistribution>,
        IMultiplyOperators<HalfCauchyDistribution, ddouble, HalfCauchyDistribution>,
        IDivisionOperators<HalfCauchyDistribution, ddouble, HalfCauchyDistribution>,
        IFittableContinuousDistribution<HalfCauchyDistribution> {

        public ddouble Gamma { get; }

        private readonly ddouble pdf_norm, gamma_inv;

        public HalfCauchyDistribution() : this(gamma: 1d) { }

        public HalfCauchyDistribution(ddouble gamma) {
            ParamAssert.ValidateScale(nameof(gamma), ParamAssert.IsFinitePositive(gamma));

            Gamma = gamma;

            pdf_norm = 2d * RcpPi / gamma;
            gamma_inv = 1d / gamma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * gamma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (IsNegative(u) || IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble pdf = pdf_norm / (Square(u) + 1d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * gamma_inv;

            if (interval == Interval.Lower) {
                if (IsNaN(u)) {
                    return NaN;
                }

                if (IsNegative(u)) {
                    return 0d;
                }
                if (IsPositiveInfinity(u)) {
                    return 1d;
                }

                ddouble cdf = 2d * RcpPi * Atan(u);

                return cdf;
            }
            else {
                if (IsNaN(u)) {
                    return NaN;
                }

                if (IsNegative(u)) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = (u < 2d)
                    ? 1d - 2d * RcpPi * Atan(u)
                    : 2d * Atan(1d / u) * RcpPi;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return 0d;
                }
                if (p >= 1d) {
                    return PositiveInfinity;
                }

                ddouble x = Gamma * TanPi(p * 0.5d);

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Gamma / TanPi(p * 0.5d);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.TanPi(u * 0.5d) * (double)Gamma;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0, PositiveInfinity);

        public override ddouble Mean => NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => Log(2d * Pi * Gamma);

        public static HalfCauchyDistribution operator *(HalfCauchyDistribution dist, ddouble k) {
            return new(dist.Gamma * k);
        }

        public static HalfCauchyDistribution operator /(HalfCauchyDistribution dist, ddouble k) {
            return new(dist.Gamma / k);
        }

        public static (HalfCauchyDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (HalfCauchyDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileScaleFitter<HalfCauchyDistribution>.Fit(new HalfCauchyDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(HalfCauchyDistribution).Name}[gamma={Gamma}]";
        }

        public override string Formula => "p(x; gamma) := 2 / (pi * (1 + u^2)) / gamma, u = x / gamma";
    }
}
