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
    public class DavisDistribution : LinearityDistribution<DavisDistribution>,
        IAdditionOperators<DavisDistribution, ddouble, DavisDistribution>,
        ISubtractionOperators<DavisDistribution, ddouble, DavisDistribution>,
        IMultiplyOperators<DavisDistribution, ddouble, DavisDistribution>,
        IDivisionOperators<DavisDistribution, ddouble, DavisDistribution>,
        IFittableContinuousDistribution<DavisDistribution> {

        public ddouble Alpha { get; }
        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, c, sigma_inv;

        private const int cache_samples = 512;
        private CDFSegmentCache? cdf_cache = null;
        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null, quantile_upper_limit_builder = null;
        private QuantileCubicApprox? sampler = null;

        public DavisDistribution(ddouble alpha) : this(alpha: alpha, mu: 0d, sigma: 1d) { }

        public DavisDistribution(ddouble alpha, ddouble sigma) : this(alpha: alpha, mu: 0d, sigma: sigma) { }

        public DavisDistribution(ddouble alpha, ddouble mu, ddouble sigma) {
            ParamAssert.ValidateShape(nameof(alpha), alpha > 1d && IsFinite(alpha));
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Alpha = alpha;
            Mu = mu;
            Sigma = sigma;

            c = 1d / (Gamma(alpha) * RiemannZeta(alpha));
            pdf_norm = c / sigma;
            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            if (u <= 0d) {
                return 0d;
            }
            if (IsNaN(u)) {
                return NaN;
            }

            ddouble pdf = pdf_norm * Pow(u, -Alpha - 1d) / Expm1(1d / u);

            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            this.cdf_cache ??= new CDFSegmentCache(
                0d, 1d,
                Integrand,
                samples: cache_samples
            );

            if (u <= 0d) {
                return interval == Interval.Lower ? 0d : 1d;
            }

            if (IsPositiveInfinity(u)) {
                return interval == Interval.Lower ? 1d : 0d;
            }

            ddouble t = 1d / (u + 1d);

            ddouble cdf = (interval == Interval.Lower)
                ? Min(1d, cdf_cache.Upper(t))
                : Min(1d, cdf_cache.Lower(t));

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

            ddouble df(ddouble u) {
                ddouble pdf = c * Pow(u, -Alpha - 1d) / Expm1(1d / u);
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return Mu;
                }

                ddouble f(ddouble u) {
                    return Min(1d, cdf_cache.Upper(1d / (u + 1d)));
                }

                this.quantile_lower_builder ??= new QuantileBuilder(
                    1d, 0d,
                    new ReadOnlyCollection<ddouble>(cdf_cache.UpperSegmentTable.Reverse().ToArray()),
                    cdf_cache.Samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_lower_builder.Estimate(p);

                if (IsNegativeInfinity(t1)) {
                    return Mu;
                }
                if (IsPositiveInfinity(t0)) {
                    return PositiveInfinity;
                }

                ddouble x = (1d - t) / t;
                ddouble x0 = (1d - t0) / t0;
                ddouble x1 = (1d - t1) / t1;

                for (int i = 0; i < 8; i++) {
                    ddouble y = f(x), dx = (y - p) / df(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x - dx, x0, x1);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                x = x * Sigma + Mu;

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                ddouble f(ddouble u) {
                    return Min(1d, cdf_cache.Lower(1d / (u + 1d)));
                }

                this.quantile_upper_builder ??= new QuantileBuilder(
                    0d, 1d,
                    cdf_cache.LowerSegmentTable,
                    cdf_cache.Samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p);

                if (IsNegativeInfinity(t0)) {
                    return PositiveInfinity;
                }
                if (IsPositiveInfinity(t1)) {
                    return Mu;
                }

                ddouble x = (1d - t) / t;
                ddouble x0 = (1d - t0) / t0;
                ddouble x1 = (1d - t1) / t1;

                if (!IsFinite(x)) {
                    ddouble r = 1d / (Alpha - 1d);

                    quantile_upper_limit_builder ??= new QuantileBuilder(
                        1024d, 0d, t => f(Pow2(t * r)), samples: 1024
                    );

                    (t, t0, t1) = quantile_upper_limit_builder.Estimate(p);

                    x = Pow2(t * r);
                    x0 = Pow2(t0 * r);
                    x1 = Pow2(t1 * r);
                }

                for (int i = 0; i < 8; i++) {
                    ddouble y = f(x), dx = (y - p) / df(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x + dx, x1, x0);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                x = x * Sigma + Mu;

                return x;
            }
        }

        internal ddouble Integrand(ddouble t) {
            if (IsZero(t)) {
                return 0d;
            }

            ddouble t_inv = 1d / t;
            ddouble u = (1d - t) * t_inv;

            ddouble pdf = c * Pow(u, -Alpha - 1d) / Expm1(1d / u);

            ddouble y = pdf * t_inv * t_inv;

            y = IsFinite(y) ? y : 0d;

            return y;
        }

        public override double Sample(Random random) {
            sampler ??= new QuantileCubicApprox(new DavisDistribution(Alpha), samples: 4096, logscale: true);

            double u = random.NextUniform();
            double r = sampler.QuantileApprox(u) * (double)Sigma + (double)Mu;

            return r;
        }

        public override (ddouble min, ddouble max) Support => (Mu, PositiveInfinity);

        public override ddouble Mean => (Alpha > 2d)
            ? Mu + Sigma * RiemannZeta(Alpha - 1d) / ((Alpha - 1d) * RiemannZeta(Alpha))
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                ddouble c = 1d / (Alpha + 1d), u = c;

                for (int i = 0; i < 256; i++) {
                    ddouble v = Exp(-1d / u);

                    ddouble du = (u * (1d - v) - c) / (1d - v * (1d + 1d / u));

                    if (!IsFinite(du)) {
                        break;
                    }

                    u = Max(u - du, u / 16d);

                    if (Abs(du) <= Abs(u) * 1e-29) {
                        break;
                    }
                }

                ddouble x = mode ??= Mu + u * Sigma;

                return x;
            }
        }

        public override ddouble Variance => (Alpha > 3d)
            ? Square(Sigma) *
                ((Alpha - 1d) * RiemannZeta(Alpha - 2d) * RiemannZeta(Alpha) - (Alpha - 2d) * Square(RiemannZeta(Alpha - 1d))) /
                ((Alpha - 2d) * Square((Alpha - 1d) * RiemannZeta(Alpha)))
            : NaN;

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            (Alpha > 4d)
            ? IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            (Alpha > 5d)
            ? IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static DavisDistribution operator +(DavisDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Sigma);
        }

        public static DavisDistribution operator -(DavisDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Sigma);
        }

        public static DavisDistribution operator *(DavisDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Sigma * k);
        }

        public static DavisDistribution operator /(DavisDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu / k, dist.Sigma / k);
        }

        public static (DavisDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (DavisDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble alpha = t / (1d - t);

                    try {
                        DavisDistribution dist = new(alpha);
                        return QuantileLinearFitter<DavisDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble alpha = t / (1d - t);
                DavisDistribution dist = new(alpha);

                return QuantileLinearFitter<DavisDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(DavisDistribution).Name}[alpha={Alpha},mu={Mu},theta={Sigma}]";
        }

        public override string Formula => "p(x; alpha, mu, sigma) := u^(- alpha - 1) / (e^(1 / u) - 1) / (gamma(alpha) * zeta(alpha) * sigma), u = (x - mu) / sigma";
    }
}
