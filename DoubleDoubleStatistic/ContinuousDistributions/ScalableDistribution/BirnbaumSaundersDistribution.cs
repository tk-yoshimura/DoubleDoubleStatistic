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
    public class BirnbaumSaundersDistribution : ScalableDistribution<BirnbaumSaundersDistribution>,
        IMultiplyOperators<BirnbaumSaundersDistribution, ddouble, BirnbaumSaundersDistribution>,
        IDivisionOperators<BirnbaumSaundersDistribution, ddouble, BirnbaumSaundersDistribution>,
        IFittableContinuousDistribution<BirnbaumSaundersDistribution> {

        public ddouble Theta { get; }
        public ddouble Alpha { get; }

        private readonly ddouble theta_inv, pdf_norm, alpha_sq;

        public BirnbaumSaundersDistribution(ddouble alpha) : this(alpha, theta: 1d) { }

        public BirnbaumSaundersDistribution(ddouble alpha, ddouble theta) {
            ParamAssert.ValidateScale(nameof(theta), ParamAssert.IsFinitePositive(theta));

            Alpha = alpha;
            Theta = theta;

            theta_inv = 1d / theta;
            alpha_sq = alpha * alpha;
            pdf_norm = theta_inv / (2d * Sqrt2 * Sqrt(PI) * alpha);
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (u <= 0d) {
                return 0d;
            }

            ddouble up1 = u + 1d, um1 = u - 1d;

            ddouble pdf = pdf_norm * up1 * Exp(-um1 * um1 / (2d * alpha_sq * u)) / ExMath.Pow3d2(u);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x)) {
                    return 1d;
                }

                ddouble cdf = Erfc((1d - u) / (Sqrt2 * Alpha * Sqrt(u))) * 0.5d;

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble cdf = Erfc((u - 1d) / (Sqrt2 * Alpha * Sqrt(u))) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p == 0d) {
                return interval == Interval.Lower ? 0 : PositiveInfinity;
            }
            if (p == 1d) {
                return interval == Interval.Lower ? PositiveInfinity : 0;
            }

            ddouble w = Alpha * InverseErfc(p * 2d), w2 = w * w;

            ddouble u;

            if (interval == Interval.Lower) {
                u = w * (w - Sqrt(w2 + 2d)) + 1d;
            }
            else {
                u = 1d / (w * (w - Sqrt(w2 + 2d)) + 1d);
            }

            ddouble x = u * Theta;

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval0();
            double w = (double)Alpha * ErrorFunction.InverseErfc(u * 2d), w2 = w * w;

            double v = (w * (w - double.Sqrt(w2 + 2d)) + 1d) * (double)Theta;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Theta * (alpha_sq + 2d) / 2d;

        public override ddouble Median => Theta;

        public override ddouble Variance => Theta * Theta * alpha_sq * (5d * alpha_sq + 4d) / 4d;

        public override ddouble Skewness => 4d * Alpha * (6d + alpha_sq * 11d) / ExMath.Pow3d2(5d * alpha_sq + 4d);

        public override ddouble Kurtosis => (48d + alpha_sq * (360d + alpha_sq * 633d)) / Square(5d * alpha_sq + 4d) - 3d;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble u = 0.25d / ExMath.Pow3d2(Alpha);

                for (int i = 0; i < 256; i++) {
                    ddouble du = (-1d + u * (-1d + 3d * alpha_sq + u * (1d + alpha_sq + u)))
                        / (-1d + 3d * alpha_sq + u * (2d * (1d + alpha_sq) + u * 3d));

                    if (!IsFinite(du)) {
                        break;
                    }

                    u -= du;

                    u = Max(Epsilon, u);

                    if (Abs(du) <= Abs(u) * 1e-29) {
                        break;
                    }
                }

                ddouble x = mode ??= u * Theta;

                return x;
            }
        }

        public static BirnbaumSaundersDistribution operator *(BirnbaumSaundersDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Theta * k);
        }

        public static BirnbaumSaundersDistribution operator /(BirnbaumSaundersDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Theta / k);
        }

        public static (BirnbaumSaundersDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BirnbaumSaundersDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble alpha = t / (1d - t);

                    try {
                        BirnbaumSaundersDistribution dist = new(alpha);
                        return QuantileScaleFitter<BirnbaumSaundersDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble alpha = t / (1d - t);
                BirnbaumSaundersDistribution dist = new(alpha);

                return QuantileScaleFitter<BirnbaumSaundersDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(BirnbaumSaundersDistribution).Name}[alpha={Alpha},theta={Theta}]";
        }

        public override string Formula => "p(x; alpha, theta) := ((u + 1) * exp(-(1 - u)^2 / (2 * alpha^2 * u))) / ((2 * u)^(3 / 2) * sqrt(pi) * alpha) / theta, u = x / theta";
    }
}
