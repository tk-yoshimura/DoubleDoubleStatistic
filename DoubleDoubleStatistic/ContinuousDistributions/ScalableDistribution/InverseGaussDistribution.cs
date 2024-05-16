using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
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
    public class InverseGaussDistribution : ScalableDistribution<InverseGaussDistribution>,
        IMultiplyOperators<InverseGaussDistribution, ddouble, InverseGaussDistribution>,
        IDivisionOperators<InverseGaussDistribution, ddouble, InverseGaussDistribution>,
        IFittableContinuousDistribution<InverseGaussDistribution> {

        public ddouble Mu { get; }
        public ddouble Lambda { get; }

        private readonly ddouble r, c, inv_mu;

        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null;

        public InverseGaussDistribution() : this(mu: 1d, lambda: 1d) { }

        public InverseGaussDistribution(ddouble mu, ddouble lambda) {
            ParamAssert.ValidateShape(nameof(mu), ParamAssert.IsFinitePositive(mu));
            ParamAssert.ValidateShape(nameof(lambda), ParamAssert.IsFinitePositive(lambda));

            Mu = mu;
            Lambda = lambda;

            c = lambda / mu;
            r = Exp(2d * c);
            inv_mu = 1d / mu;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * inv_mu;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble pdf = Sqrt(c / (2d * PI * u)) * Exp(c * Square(u - 1d) / (-2d * u)) / u * inv_mu;
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * inv_mu;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble v = Sqrt(c / (u * 2d));
            ddouble vn = v * (1d - u), vp = v * (1d + u);

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }
                if (IsPositiveInfinity(u)) {
                    return 1d;
                }

                ddouble cdf = (Erfc(vn) + r * Erfc(vp)) * 0.5d;
                cdf = Min(cdf, 1d);

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = (Erfc(-vn) - r * Erfc(vp)) * 0.5d;
                cdf = Max(cdf, 0d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            ddouble df(ddouble x) {
                ddouble y = Sqrt(c / (2d * PI * x)) * Exp(c * Square(x - 1d) / (-2d * x)) / x;
                return y;
            }

            if (interval == Interval.Lower) {
                ddouble f(ddouble x) {
                    ddouble u = Sqrt(c / (x * 2d));
                    ddouble y = (r * Erfc(u * (x + 1d)) + Erfc(u * (1d - x))) * 0.5d;
                    return Min(1d, y);
                }

                this.quantile_lower_builder ??= new QuantileBuilder(0d, 2d, f, samples: 1024);

                (ddouble x, ddouble x0, ddouble x1) = quantile_lower_builder.Estimate(p);

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

                x *= Mu;

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                ddouble f(ddouble x) {
                    ddouble u = Sqrt(c / (x * 2d));
                    ddouble y = (-r * Erfc(u * (x + 1d)) + Erfc(u * (x - 1d))) * 0.5d;
                    return Max(0d, y);
                }

                this.quantile_upper_builder ??= new QuantileBuilder(0d, 1d,
                    t => t > 0d ? f((1d - t) / t) : 0d,
                    samples: 1024
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p);

                ddouble x = (1d - t) / t;
                ddouble x0 = (1d - t0) / t0;
                ddouble x1 = (1d - t1) / t1;

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

                x *= Mu;

                return x;
            }
        }

        public override double Sample(Random random) {
            double lambda = (double)Lambda, mu = (double)Mu;

            double x = random.NextGaussian();
            double y = x * x * mu;
            double z = random.NextUniform();
            double w = mu - (0.5 * mu / lambda) * (double.Sqrt(y * (y + 4d * lambda)) - y);

            double r = (z < (mu / (mu + w))) ? w : (mu * mu / w);

            return r;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Mu;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode {
            get {
                ddouble c = (3d * Mu) / (2d * Lambda);
                return Mu * (Sqrt(1d + Square(c)) - c);
            }
        }

        public override ddouble Variance => Cube(Mu) / Lambda;

        public override ddouble Skewness => 3d * Sqrt(Mu / Lambda);

        public override ddouble Kurtosis => 15d * Mu / Lambda;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static InverseGaussDistribution operator *(InverseGaussDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Lambda * k);
        }

        public static InverseGaussDistribution operator /(InverseGaussDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Lambda / k);
        }

        public static (InverseGaussDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (InverseGaussDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble mu = t / (1d - t);

                    try {
                        InverseGaussDistribution dist = new(mu, 1d);
                        return QuantileScaleFitter<InverseGaussDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble mu = t / (1d - t);
                InverseGaussDistribution dist = new(mu, 1d);

                return QuantileScaleFitter<InverseGaussDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(InverseGaussDistribution).Name}[mu={Mu},lambda={Lambda}]";
        }

        public override string Formula => "p(x; mu, lambda) := sqrt(c / (2 * pi * u^3)) * exp(-c * (u - 1)^2 / (2 * u)) / mu, c = lambda / mu, u = x / mu";
    }
}
