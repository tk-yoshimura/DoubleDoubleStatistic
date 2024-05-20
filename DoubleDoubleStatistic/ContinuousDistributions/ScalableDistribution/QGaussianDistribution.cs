using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class QGaussianDistribution : ScalableDistribution<QGaussianDistribution>,
        IMultiplyOperators<QGaussianDistribution, ddouble, QGaussianDistribution>,
        IDivisionOperators<QGaussianDistribution, ddouble, QGaussianDistribution>,
        IFittableContinuousDistribution<QGaussianDistribution> {

        public ddouble Q { get; }
        public ddouble Sigma { get; }

        private readonly ddouble cq, power, sigma_inv;
        private CDFSegmentCache? cdf_cache;
        private QuantileBuilder? quantile_lower_builder = null;

        public QGaussianDistribution(ddouble q) : this(q, sigma: 1d) { }

        public QGaussianDistribution(ddouble q, ddouble sigma) {
            ParamAssert.ValidateShape(nameof(q), q >= 0 && q < 3);
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Q = q;
            Sigma = sigma;

            sigma_inv = 1d / sigma;

            ddouble sqrt_pi = Sqrt(PI);

            if (q < 1d) {
                cq = 2d * sqrt_pi * Gamma(1d / (1d - q)) /
                    ((3d - q) * Sqrt(1d - q) * Gamma((3d - q) / (2d * (1d - q))));
            }
            else if (q == 1d) {
                cq = sqrt_pi;
            }
            else {
                cq = sqrt_pi * Gamma((3d - q) / (2d * (q - 1d))) /
                    (Sqrt(q - 1d) * Gamma(1d / (q - 1d)));
            }

            power = 1d / (1d - q);
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (IsInfinity(u)) {
                return 0d;
            }

            if (IsFinite(power)) {
                ddouble pdf = sigma_inv / cq * Pow(Max(0d, 1d + (Q - 1d) * u * u), power);

                if (IsNaN(pdf)) {
                    return 0d;
                }

                return pdf;
            }
            else {
                ddouble pdf = sigma_inv / cq * Exp(-u * u);

                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsPositive(x) && interval == Interval.Lower) {
                return 1d - CDF(x, Interval.Upper);
            }
            if (IsNegative(x) && interval == Interval.Upper) {
                return 1d - CDF(x, Interval.Lower);
            }

            ddouble u = Abs(x) * sigma_inv;

            if (Abs(u) < 1e-30) {
                return 0.5d;
            }

            if (IsNaN(u)) {
                return NaN;
            }

            if (!IsFinite(power)) {
                ddouble cdf = Erfc(u);

                return cdf * 0.5d;
            }
            else {
                cdf_cache ??= SetupCDFCache();

                if (Q < 1) {
                    ddouble range = 1d / Sqrt(1d - Q);

                    if (u > range) {
                        return IsNegative(x) ^ interval == Interval.Lower ? 1d : 0d;
                    }

                    ddouble cdf = cdf_cache.Upper(u);

                    return cdf;
                }
                else {
                    ddouble t = 1d / (Abs(u) + 1d);

                    ddouble cdf = cdf_cache.Lower(t);

                    return cdf;
                }
            }
        }

        private CDFSegmentCache SetupCDFCache() {
            if (Q < 1) {
                ddouble range = 1d / Sqrt(1d - Q);

                return new CDFSegmentCache(
                    0d, range,
                    u => Pow(Max(0d, 1d + (Q - 1d) * u * u), power) / cq,
                    samples: 4096
                );
            }
            else {
                return new CDFSegmentCache(
                    0d, 1d,
                    t => {
                        ddouble t_inv = 1d / t;
                        ddouble u = (1d - t) * t_inv;

                        ddouble pdf = Pow(Max(0d, 1d + (Q - 1d) * u * u), power) / cq;

                        ddouble y = pdf * t_inv * t_inv;

                        y = IsFinite(y) ? y : 0d;

                        return y;
                    },
                    samples: 16384
                );
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (!IsFinite(power)) {
                ddouble x = Sigma * InverseErfc(2d * p);

                x = (interval == Interval.Lower) ? -x : x;

                return x;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            if (interval == Interval.Upper) {
                return -Quantile(p, Interval.Lower);
            }

            if (Abs(p - 0.5d) < 1e-30) {
                return 0d;
            }

            cdf_cache ??= SetupCDFCache();

            if (Q < 1) {
                ddouble range = 1d / Sqrt(1d - Q);

                if (p <= 0d) {
                    return -Sigma * range;
                }

                this.quantile_lower_builder ??= new QuantileBuilder(
                    0d, range,
                    new ReadOnlyCollection<ddouble>(cdf_cache.LowerSegmentTable.Select(x => Min(0.5d, x)).ToArray()),
                    cdf_cache.Samples
                );

                (ddouble x, ddouble x0, ddouble x1) = quantile_lower_builder.Estimate(0.5d - p);

                x *= -Sigma;
                x0 *= -Sigma;
                x1 *= -Sigma;

                for (int i = 0; i < 64; i++) {
                    ddouble y = CDF(x, Interval.Lower), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x - dx, x1, x0);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }

            else {
                if (p <= 0d) {
                    return NegativeInfinity;
                }

                this.quantile_lower_builder ??= new QuantileBuilder(
                    1d, 0d,
                    new ReadOnlyCollection<ddouble>(cdf_cache.UpperSegmentTable.Reverse().Select(x => Min(0.5d, x)).ToArray()),
                    cdf_cache.Samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_lower_builder.Estimate(0.5d - p);

                if (IsNegativeInfinity(t1)) {
                    return NegativeInfinity;
                }
                if (IsPositiveInfinity(t0)) {
                    return 0d;
                }

                ddouble x = (1d - t) / -t * Sigma;
                ddouble x0 = (1d - t0) / -t0 * Sigma;
                ddouble x1 = (1d - t1) / -t1 * Sigma;

                for (int i = 0; i < 64; i++) {
                    ddouble y = CDF(x, Interval.Lower), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x - dx, x1, x0);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }
        }

        public override double Sample(Random random) {
            if (Q == 1d) {
                double w = random.NextGaussian() * double.Sqrt(0.5d) * (double)Sigma;

                return w;
            }
            else {
                double q = (double)Q, q_prime = 1d - (1d + q) / (3d - q);
                double u = random.NextUniformOpenInterval0(), theta = random.NextUniform();

                double r = double.Sqrt(-2d * (double.Pow(u, q_prime) - 1) / q_prime) * double.CosPi(2 * theta) / double.Sqrt(3d - q);

                double w = r * (double)Sigma;

                return w;
            }
        }

        public override (ddouble min, ddouble max) Support => Q < 1d
            ? (-Sigma / Sqrt(1d - Q), Sigma / Sqrt(1d - Q))
            : (NegativeInfinity, PositiveInfinity);

        public override ddouble Mean => (Q < 2d)
            ? 0d
            : NaN;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Variance => (Q * 3d < 5d)
            ? Sigma * Sigma / (5d - 3d * Q)
            : (Q < 2d) ? PositiveInfinity : NaN;

        public override ddouble Skewness => (Q < 1.5d)
            ? 0d
            : NaN;

        public override ddouble Kurtosis => (Q * 5d < 7d)
            ? 6d * (Q - 1d) / (7d - 5d * Q)
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 32768);

        public static QGaussianDistribution operator *(QGaussianDistribution dist, ddouble k) {
            return new(dist.Q, dist.Sigma * k);
        }

        public static QGaussianDistribution operator /(QGaussianDistribution dist, ddouble k) {
            return new(dist.Q, dist.Sigma / k);
        }

        public static (QGaussianDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (QGaussianDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble q = BisectionMinimizeSearch1D.Search(
                q => {
                    try {
                        QGaussianDistribution dist = new(q, 1d);
                        return QuantileScaleFitter<QGaussianDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (0, 2.9995d), iter: 32
            );

            try {
                QGaussianDistribution dist = new(q, 1d);

                return QuantileScaleFitter<QGaussianDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(QGaussianDistribution).Name}[q={Q},sigma={Sigma}]";
        }

        public override string Formula => "p(x; q, sigma) := e_q(-u^2) / (sigma * c_q), u = x / sigma";
    }
}
