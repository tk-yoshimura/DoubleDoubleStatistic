using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class AlphaDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<AlphaDistribution> {

        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm, cdf_norm;

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        public AlphaDistribution() : this(alpha: 1d) { }

        public AlphaDistribution(ddouble alpha) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));

            Alpha = alpha;
            cdf_norm = 1d / Erfc(-alpha * sqrt2_inv);
            pdf_norm = cdf_norm * Sqrt(2d * RcpPI);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble u = Alpha - 1d / x;

            ddouble pdf = pdf_norm / (x * x) * Exp(u * u * -0.5d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = Alpha - 1d / x;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x)) {
                    return 1d;
                }

                ddouble cdf = cdf_norm * Erfc(-u * sqrt2_inv);

                if (IsNaN(cdf)) {
                    return 1d;
                }

                cdf = Min(1d, cdf);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble cdf = 1d - cdf_norm * Erfc(-u * sqrt2_inv);

                if (IsNaN(cdf)) {
                    return 0d;
                }

                cdf = Max(0d, cdf);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p, Interval.Lower);
            }

            if (p == 1d) {
                return PositiveInfinity;
            }

            ddouble u = InverseErfc(p / cdf_norm) * Sqrt2;
            ddouble x = 1d / Max(0d, Alpha + u);

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = ErrorFunction.InverseErfc(u / (double)cdf_norm) * double.Sqrt(2d);
            double w = 1d / double.Max(0d, (double)Alpha + v);

            return w;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => (Sqrt(Alpha * Alpha + 8d) - Alpha) * 0.25d;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (AlphaDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (AlphaDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble alpha = t / (1d - t);

                    try {
                        AlphaDistribution dist = new(alpha);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble alpha = t / (1d - t);
                AlphaDistribution dist = new(alpha);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(AlphaDistribution).Name}[alpha={Alpha}]";
        }

        public override string Formula => "p(x; alpha) := sqrt(2 / pi) / erfc(-alpha / sqrt(2)) * exp(-(alpha - 1 / x)^2 / 2) / x^2";
    }
}
