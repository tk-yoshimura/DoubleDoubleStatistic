using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class TukeyLambdaDistribution : ContinuousDistribution { //,
        //IFittableContinuousDistribution<TukeyLambdaDistribution> {

        public ddouble Lambda { get; }

        private readonly ReadOnlyCollection<ddouble> quantile_table;
        private const int quantile_samples = 256;

        public TukeyLambdaDistribution(ddouble lambda) {
            ParamAssert.ValidateShape(nameof(lambda), IsFinite(lambda));

            Lambda = lambda;

            ddouble[] quantile_table = new ddouble[quantile_samples + 1];
            quantile_table[0] = lambda > 0d ? 1d / lambda : PositiveInfinity;
            for (int i = 1; i < quantile_samples; i++) {
                quantile_table[i] = -Quantile((ddouble)i / (quantile_samples * 2), Interval.Lower);
            }
            quantile_table[^1] = 0d;

            this.quantile_table = new ReadOnlyCollection<ddouble>(quantile_table.Reverse().ToArray());
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            x = Abs(x);

            if ((Lambda > 0d && x * Lambda > 1d) || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble cdf = CDF(x, Interval.Upper);

            ddouble pdf = 1d / (Pow(cdf, Lambda - 1d) + Pow(1d - cdf, Lambda - 1d));

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) && interval == Interval.Upper) {
                return 1d - CDF(-x, Interval.Upper);
            }
            if (IsPositive(x) && interval == Interval.Lower) {
                return 1d - CDF(-x, Interval.Lower);
            }

            x = Abs(x);

            if ((Lambda > 0d && x * Lambda > 1d) || IsPositiveInfinity(x)) {
                return 0d;
            }

            int index = Indexer.BisectionSearch(x, quantile_table);

            ddouble x0 = quantile_table[index], x1 = quantile_table[index + 1];
            ddouble p0 = (ddouble)(quantile_samples - index - 1) / (quantile_samples * 2);
            ddouble p1 = (ddouble)(quantile_samples - index) / (quantile_samples * 2);

            ddouble p = (x - x0) / (x1 - x0) * (p1 - p0) + p0;

            for (int i = 0; i < 8; i++) {
                ddouble q = (Lambda != 0d) 
                    ? (Pow(p, Lambda) - Pow(1d - p, Lambda)) / Lambda
                    : Log(p / (1d - p));

                ddouble dq = Pow(p, Lambda - 1d) + Pow(1d - p, Lambda - 1d);

                ddouble dp = (x + q) / dq;

                if (!IsFinite(dp)) {
                    break;
                }

                p = Clamp(p - dp, p0, p1);

                if (Abs(dp) <= Abs(p) * 1e-29) {
                    break;
                }
            }

            return p;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return -Quantile(p, Interval.Lower);
            }

            if (Lambda != 0d) {
                ddouble x = (Pow(p, Lambda) - Pow(1d - p, Lambda)) / Lambda;

                return x;
            }
            else {
                return Log(p / (1d - p));
            }
        }

        public override double Sample(Random random) {
            double p = random.NextUniformOpenInterval01();
            double lambda = (double)Lambda;

            if (lambda != 0d) {
                double x = (double.Pow(p, lambda) - double.Pow(1d - p, lambda)) / lambda;

                return x;
            }
            else {
                return double.Log(p / (1d - p));
            }
        }

        public override (ddouble min, ddouble max) Support => Lambda > 0d
            ? (-1d / Lambda, 1d / Lambda)
            : (NegativeInfinity, PositiveInfinity);

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Mean => (Lambda > -1d) ? 0d : NaN;

        public override ddouble Variance => (Lambda > -0.5d)
            ? 2d * (1d / (1d + 2d * Lambda) - Square(Gamma(Lambda + 1d)) / Gamma(2d * Lambda + 2d)) / Square(Lambda)
            : Square(PI) / 3d;

        public override ddouble Skewness => (Lambda * 3d > -1d) ? 0d : NaN;

        public override ddouble Kurtosis {
            get {
                if (Lambda * 4d < -1d) {
                    return NaN;
                }

                ddouble g1 = Gamma(Lambda + 1d), g2 = Gamma(2d * Lambda + 1d);
                ddouble g3 = Gamma(3d * Lambda + 1d), g4 = Gamma(4d * Lambda + 1d);

                if (Lambda == 0d) {
                    return (ddouble)6 / 5;
                }

                ddouble kurtosis = Square(g2 * (2d * Lambda + 1d) / (g1 * g1 - g2)) * (3d * g2 * g2 - 4d * g1 * g3 + g4) 
                    / ((8d * Lambda + 2d) * g4) - 3d;

                return kurtosis;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (TukeyLambdaDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (TukeyLambdaDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble lambda = BisectionMinimizeSearch1D.Search(
                lambda => {
                    try {
                        TukeyLambdaDistribution dist = new(lambda);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (-8d, 8d), iter: 32
            );

            try {
                TukeyLambdaDistribution dist = new(lambda);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(RiceDistribution).Name}[lambda={Lambda}]";
        }

        public override string Formula => "p(x; lambda)";
    }
}
