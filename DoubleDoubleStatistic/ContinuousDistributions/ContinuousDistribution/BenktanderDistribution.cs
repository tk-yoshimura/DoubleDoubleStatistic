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
    public class BenktanderDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<BenktanderDistribution> {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble c;

        private QuantileBuilder? quantile_upper_builder = null;

        private QuantileSampler? sampler = null;

        public BenktanderDistribution(ddouble alpha, ddouble beta) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateShape(nameof(beta), beta > 0d && beta <= alpha * (alpha + 1d) * 0.5d);

            Alpha = alpha;
            Beta = beta;

            c = beta / alpha;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x < 1d || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble lnx = Log(x), beta_lnx = Beta * lnx;

            ddouble pdf = Exp(lnx * (-Alpha - 2d - beta_lnx)) * ((Alpha + 2d * beta_lnx + 1d) * (2d * c * lnx + 1d) - 2d * c);

            pdf = Max(0d, pdf);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return 1d - CDF(x, Interval.Upper);
            }
            else {
                if (x <= 1d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble lnx = Log(x), beta_lnx = Beta * lnx;

                ddouble cdf = Exp(lnx * (-Alpha - 1d - beta_lnx)) * (2d * c * lnx + 1d);

                cdf = Clamp(cdf, 0d, 1d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return Quantile(1d - p, Interval.Upper);
            }

            if (p <= 0) {
                return PositiveInfinity;
            }
            if (p >= 1) {
                return 0d;
            }

            this.quantile_upper_builder ??= new QuantileBuilder(0d, 1d,
                t => t > 0d ? CDF((1d - t) / t, Interval.Upper) : 0d,
                samples: 1024
            );

            (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p);

            ddouble x = (1d - t) / t;
            ddouble x0 = (1d - t0) / t0;
            ddouble x1 = (1d - t1) / t1;

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

            x = Max(1d, x);

            return x;
        }

        public override double Sample(Random random) {
            sampler ??= new QuantileSampler(this, samples: 4096);

            double u = random.NextUniform();
            double r = sampler.QuantileApprox(u);

            return r;
        }

        public override (ddouble min, ddouble max) Support => (1d, PositiveInfinity);

        public override ddouble Mean => 1d + 1d / Alpha;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode =>
            Max(1d, Exp((-Alpha - 1d + Sqrt(6d * Beta + 1d)) / (2d * Beta)));

        public override ddouble Variance =>
            Sqrt(PI) * Exp(Square(Alpha - 1d) / (4d * Beta)) * Erfc((Alpha - 1d) / (2d * Sqrt(Beta))) / (Alpha * Sqrt(Beta)) - 1d / Square(Alpha);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (BenktanderDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BenktanderDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble mean = samples.Mean();
            ddouble alpha = Max(1e-2, 1d / (mean - 1d));

            ddouble beta_max = alpha * (alpha + 1d) * 0.5d;

            ddouble beta = BisectionMinimizeSearch1D.Search(
                beta => {
                    try {
                        BenktanderDistribution dist = new(alpha, beta);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-4 * beta_max, 0.9999d * beta_max), iter: 32
            );

            try {
                BenktanderDistribution dist = new(alpha, beta);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(BenktanderDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := x^(- alpha - 2) * exp(-beta * log(x)^2) * ((alpha + 2 * beta * log(x) + 1) * (2 * beta * log(x) / alpha + 1) - 2 * beta / alpha)";
    }
}
