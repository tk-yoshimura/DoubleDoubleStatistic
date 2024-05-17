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
    public class DagumDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<DagumDistribution> {

        public ddouble A { get; }
        public ddouble P { get; }

        private readonly ddouble ap, a_inv, p_inv;

        public DagumDistribution(ddouble a, ddouble p) {
            ParamAssert.ValidateShape(nameof(a), ParamAssert.IsFinitePositive(a));
            ParamAssert.ValidateShape(nameof(p), ParamAssert.IsFinitePositive(p));

            A = a;
            P = p;

            ap = a * p;
            a_inv = 1d / a;
            p_inv = 1d / p;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (A == 1d) {
                ddouble pdf = P * Pow(x, P - 1d) / Pow(x + 1d, P + 1d);
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }
            else {
                ddouble xa = Pow(x, A);

                if (xa <= 0d) {
                    return A < 1d ? PositiveInfinity : 0d;
                }

                if (IsPositiveInfinity(xa)) {
                    return 0d;
                }

                ddouble pdf = ap * Pow(xa, P) / (x * Pow(xa + 1d, P + 1d));
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble xa = Pow(x, -A), xcp1 = 1d + xa;

            if (interval == Interval.Lower) {
                if (IsNegative(x) || xa <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 1d;
                }

                ddouble cdf = Min(1d, Pow(xcp1, -P));

                return cdf;
            }
            else {
                if (IsNegative(x) || xa <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 0d;
                }

                ddouble cdf = Max(0d, 1d - Pow(xcp1, -P));

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

                ddouble x = Pow(Pow(1d / p, p_inv) - 1d, -a_inv);

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Pow(Pow(1d / (1d - p), p_inv) - 1d, -a_inv);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(double.Pow(1d / u, (double)p_inv) - 1d, -(double)a_inv);

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            A > 1d
            ? P * Beta(P + a_inv, 1d - a_inv)
            : NaN;


        public override ddouble Median =>
            Pow(Pow2(p_inv) - 1d, -a_inv);

        public override ddouble Mode =>
            P * A > 1d
            ? Pow((A + 1d) / (P * A - 1d), -a_inv)
            : 0d;

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            A > 2d
            ? IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            A > 3d
            ? IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            A > 4d
            ? IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (DagumDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (DagumDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble a = t.u / (1d - t.u);
                    ddouble p = t.v / (1d - t.v);

                    try {
                        DagumDistribution dist = new(a, p);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-10d, 1e-10d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble a = u / (1d - u);
                ddouble p = v / (1d - v);
                DagumDistribution dist = new(a, p);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(DagumDistribution).Name}[a={A},p={P}]";
        }

        public override string Formula => "p(x; a, p) := a * p * x^(a * p - 1) / (1 + x^a)^(p + 1)";
    }
}
