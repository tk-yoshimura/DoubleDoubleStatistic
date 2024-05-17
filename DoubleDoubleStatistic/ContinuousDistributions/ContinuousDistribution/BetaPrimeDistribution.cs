using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BetaPrimeDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<BetaPrimeDistribution> {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble pdf_lognorm;

        private readonly BetaDistribution randam_gen_beta_dist;

        public BetaPrimeDistribution(ddouble alpha, ddouble beta) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateShape(nameof(beta), ParamAssert.IsFinitePositive(beta));

            Alpha = alpha;
            Beta = beta;

            pdf_lognorm = LogBeta(alpha, beta) * LbE;

            randam_gen_beta_dist = new(alpha, beta);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || IsPositiveInfinity(x)) {
                return 0d;
            }

            if (x <= 0d) {
                return Alpha < 1d ? PositiveInfinity : Alpha == 1d ? Beta : 0d;
            }

            ddouble pdf = Pow2(Log2(x) * (Alpha - 1d) - Log1p(x) * (Alpha + Beta) * LbE - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble xp1 = x + 1d;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(xp1)) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(x / xp1, Alpha, Beta);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(1d / xp1, Beta, Alpha);

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

                ddouble u = InverseIncompleteBeta(p, Alpha, Beta);
                ddouble x = u / (1d - u);

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble u = InverseIncompleteBeta(1d - p, Alpha, Beta);
                ddouble x = u / (1d - u);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = randam_gen_beta_dist.Sample(random);

            double r = u / (1d - u);

            return r;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Beta > 1d ? Alpha / (Beta - 1d) : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Alpha >= 1d ?
            (Alpha - 1d) / (Beta + 1d)
            : 0d;

        public override ddouble Variance => Beta > 2d
            ? Alpha * (Alpha + Beta - 1d) / ((Beta - 2d) * Square(Beta - 1d))
            : NaN;

        public override ddouble Skewness => Beta > 3d
            ? 2d * (2d * Alpha + Beta - 1d) / (Beta - 3d) * Sqrt((Beta - 2d) / (Alpha * (Alpha + Beta - 1d)))
            : NaN;

        public override ddouble Kurtosis => Beta > 4d
            ? 6d * (Alpha * (Alpha + Beta - 1d) * (5d * Beta - 11d) + Square(Beta - 1d) * (Beta - 2d)) / (Alpha * (Alpha + Beta - 1d) * (Beta - 3d) * (Beta - 4d))
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (BetaPrimeDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BetaPrimeDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble alpha = t.u / (1d - t.u);
                    ddouble beta = t.v / (1d - t.v);

                    try {
                        BetaPrimeDistribution dist = new(alpha, beta);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-10d, 1e-10d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble alpha = u / (1d - u);
                ddouble beta = v / (1d - v);
                BetaPrimeDistribution dist = new(alpha, beta);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(BetaPrimeDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := x^(alpha - 1) * (1 + x)^(- alpha - beta) / beta(alpha, beta)";
    }
}
