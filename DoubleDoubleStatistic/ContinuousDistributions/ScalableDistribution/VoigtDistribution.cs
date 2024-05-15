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
using Complex = DoubleDoubleComplex.Complex;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class VoigtDistribution : ScalableDistribution<VoigtDistribution>,
        IAdditionOperators<VoigtDistribution, VoigtDistribution, VoigtDistribution>,
        IMultiplyOperators<VoigtDistribution, ddouble, VoigtDistribution>,
        IDivisionOperators<VoigtDistribution, ddouble, VoigtDistribution>,
        IFittableContinuousDistribution<VoigtDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private readonly ddouble pdf_norm, c, sigma_inv, zr, cdf_limit;

        public VoigtDistribution() : this(gamma: 1d, sigma: 1d) { }

        private const int cache_samples = 512;
        private CDFSegmentCache? cdf_cache;
        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null;

        public VoigtDistribution(ddouble gamma, ddouble sigma) {
            ValidateScale(gamma);
            ValidateScale(sigma);

            Gamma = gamma;
            Sigma = sigma;

            pdf_norm = 1d / (sigma * Sqrt(2d * PI));
            c = 1d / Sqrt(2d * PI);
            sigma_inv = 1d / sigma;
            zr = gamma / (Sqrt2 * sigma);
            cdf_limit = gamma / sigma * RcpPI;
        }

        public override ddouble PDF(ddouble x) {
            if (IsInfinity(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = x * sigma_inv;

            Complex z = (zr, -u * sqrt2_inv);

            ddouble pdf = Complex.Erfcx(z).R * pdf_norm;
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsInfinity(x)) {
                if (x > 0d) {
                    return interval == Interval.Lower ? 1d : 0d;
                }
                else {
                    return interval == Interval.Lower ? 0d : 1d;
                }
            }

            ddouble p = PDF(x);
            if (IsZero(p)) {
                return x < 0d ? 0d : 1d;
            }

            this.cdf_cache ??= new CDFSegmentCache(
                0d, 1d,
                Integrand,
                samples: cache_samples
            );

            ddouble u = x * sigma_inv;

            ddouble t = 1d / (Abs(u) + 1d);

            ddouble cdf = (interval == Interval.Lower) ^ (u < 0d)
                ? Min(1d, 0.5d + cdf_cache.Upper(t))
                : Min(0.5d, cdf_cache.Lower(t));

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            this.cdf_cache ??= new CDFSegmentCache(
                0d, 1d,
                Integrand,
                samples: cache_samples
            );

            if (interval == Interval.Lower) {
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

                for (int i = 0; i < 8; i++) {
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
                    return PositiveInfinity;
                }

                this.quantile_upper_builder ??= new QuantileBuilder(
                    0d, 1d,
                    new ReadOnlyCollection<ddouble>(cdf_cache.LowerSegmentTable.Select(x => Min(0.5d, x)).ToArray()),
                    cdf_cache.Samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p);

                if (IsNegativeInfinity(t0)) {
                    return PositiveInfinity;
                }
                if (IsPositiveInfinity(t1)) {
                    return 0d;
                }

                ddouble x = (1d - t) / t * Sigma;
                ddouble x0 = (1d - t0) / t0 * Sigma;
                ddouble x1 = (1d - t1) / t1 * Sigma;

                for (int i = 0; i < 8; i++) {
                    ddouble y = CDF(x, Interval.Upper), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x + dx, x1, x0);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01() - 0.5d;
            double r = double.TanPi(u) * (double)Gamma;
            double v = random.NextGaussian() * (double)Sigma;

            double w = r + v;

            return w;
        }

        internal ddouble Integrand(ddouble t) {
            if (IsZero(t)) {
                return cdf_limit;
            }

            ddouble t_inv = 1d / t;
            ddouble u = (1d - t) * t_inv;

            Complex z = (zr, -u * sqrt2_inv);
            ddouble pdf = Complex.Erfcx(z).R * c;

            ddouble y = pdf * t_inv * t_inv;

            y = IsFinite(y) ? y : 0d;

            return y;
        }

        public override bool AdditiveClosed => true;

        public override bool Symmetric => true;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static VoigtDistribution operator +(VoigtDistribution dist1, VoigtDistribution dist2) {
            return new(dist1.Gamma + dist2.Gamma, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static VoigtDistribution operator *(VoigtDistribution dist, ddouble k) {
            return new(dist.Gamma * k, dist.Sigma * k);
        }

        public static VoigtDistribution operator /(VoigtDistribution dist, ddouble k) {
            return new(dist.Gamma / k, dist.Sigma / k);
        }

        public static (VoigtDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (VoigtDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble gamma = t / (1d - t);

                    try {
                        VoigtDistribution dist = new(gamma, 1d);
                        return QuantileScaleFitter<VoigtDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10, 1000d / 1001d), iter: 32
            );

            try {
                ddouble gamma = t / (1d - t);
                VoigtDistribution dist = new(gamma, 1d);

                return QuantileScaleFitter<VoigtDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string Formula => "p(x; gamma, sigma) := Re[erfcx((x + i * gamma) / (sigma * sqrt(2)))] / (sigma * sqrt(2 * pi))";

        public override string ToString() {
            return $"{typeof(VoigtDistribution).Name}[gamma={Gamma},sigma={Sigma}]";
        }
    }
}
