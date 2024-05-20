using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class SkewCauchyDistribution : LinearityDistribution<SkewCauchyDistribution>,
        IAdditionOperators<SkewCauchyDistribution, ddouble, SkewCauchyDistribution>,
        ISubtractionOperators<SkewCauchyDistribution, ddouble, SkewCauchyDistribution>,
        IMultiplyOperators<SkewCauchyDistribution, ddouble, SkewCauchyDistribution>,
        IDivisionOperators<SkewCauchyDistribution, ddouble, SkewCauchyDistribution>,
        IFittableContinuousDistribution<SkewCauchyDistribution> {

        public ddouble Mu { get; }
        public ddouble Gamma { get; }
        public ddouble Alpha { get; }

        private readonly ddouble gamma_inv, alphap1, alpham1, u0p, alphav;

        public SkewCauchyDistribution(ddouble alpha) : this(alpha, mu: 0d, gamma: 1d) { }
        public SkewCauchyDistribution(ddouble alpha, ddouble sigma) : this(alpha, mu: 0d, gamma: sigma) { }

        public SkewCauchyDistribution(ddouble alpha, ddouble mu, ddouble gamma) {
            ParamAssert.ValidateShape(nameof(alpha), Abs(alpha) <= 1d);
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(gamma), ParamAssert.IsFinitePositive(gamma));

            Mu = mu;
            Gamma = gamma;
            Alpha = alpha;

            gamma_inv = 1d / gamma;
            alphap1 = 1d + Alpha;
            alpham1 = 1d - Alpha;

            u0p = CDF(Mu);

            alphav = Sqrt(1d / (1d + alpha * alpha));
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * gamma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble pdf = 1d / (PI * (u * u / Square(Alpha * Sign(u) + 1d) + 1d));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * gamma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (IsInfinity(u)) {
                return (IsNegative(u) ^ interval == Interval.Lower) ? 1d : 0d;
            }

            if (interval == Interval.Lower) {
                if (u < 0d) {
                    ddouble cdf = alpham1 * 0.5d + alpham1 * Atan(u / alpham1) * RcpPI;

                    return cdf;
                }
                else {
                    ddouble cdf = alpham1 * 0.5d + alphap1 * Atan(u / alphap1) * RcpPI;

                    return cdf;
                }
            }
            else {
                if (u < 0d) {
                    ddouble cdf = alphap1 * 0.5d - alpham1 * Atan(u / alpham1) * RcpPI;

                    return cdf;
                }
                else {
                    ddouble cdf = alphap1 * 0.5d - alphap1 * Atan(u / alphap1) * RcpPI;

                    return cdf;
                }
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble u;

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return NegativeInfinity;
                }
                if (p >= 1d) {
                    return PositiveInfinity;
                }

                if (p < u0p) {
                    u = TanPI((p - alpham1 * 0.5d) / alpham1) * alpham1;
                }
                else {
                    u = TanPI((p - alpham1 * 0.5d) / alphap1) * alphap1;
                }
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return NegativeInfinity;
                }

                if (1d - p < u0p) {
                    u = -TanPI((p - alphap1 * 0.5d) / alpham1) * alpham1;
                }
                else {
                    u = -TanPI((p - alphap1 * 0.5d) / alphap1) * alphap1;
                }
            }

            ddouble x = u * Gamma + Mu;

            return x;
        }

        public override double Sample(Random random) {
            double p = random.NextUniformOpenInterval01();

            double u;
            if (p < (double)u0p) {
                u = double.TanPi((p - (double)alpham1 * 0.5d) / (double)alpham1) * (double)alpham1;
            }
            else {
                u = double.TanPi((p - (double)alpham1 * 0.5d) / (double)alphap1) * (double)alphap1;
            }

            double w = u * (double)Gamma + (double)Mu;

            return w;
        }

        public override ddouble Mean => NaN;

        public override ddouble Mode => throw new NotFiniteNumberException();

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 32768);

        public static SkewCauchyDistribution operator +(SkewCauchyDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Gamma);
        }

        public static SkewCauchyDistribution operator -(SkewCauchyDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Gamma);
        }

        public static SkewCauchyDistribution operator *(SkewCauchyDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Gamma * k);
        }

        public static SkewCauchyDistribution operator /(SkewCauchyDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu / k, dist.Gamma / k);
        }

        public static (SkewCauchyDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (SkewCauchyDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble alpha = BisectionMinimizeSearch1D.Search(
                alpha => {
                    try {
                        SkewCauchyDistribution dist = new(alpha);
                        return QuantileLinearFitter<SkewCauchyDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (-1d, 1d), iter: 32
            );

            try {
                SkewCauchyDistribution dist = new(alpha);

                return QuantileLinearFitter<SkewCauchyDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(SkewCauchyDistribution).Name}[alpha={Alpha},mu={Mu},gamma={Gamma}]";
        }

        public override string Formula => "p(x; alpha, mu, gamma) := 1 / (pi * (u^2 / (alpha * sign(x) + 1)^2 + 1)), u = (x - mu) / gamma";
    }
}
