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
    public class KumaraswamyDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<KumaraswamyDistribution> {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble ab, alpha_inv, beta_inv;

        public KumaraswamyDistribution(ddouble alpha, ddouble beta) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateShape(nameof(beta), ParamAssert.IsFinitePositive(beta));

            Alpha = alpha;
            Beta = beta;

            ab = alpha * beta;
            alpha_inv = 1d / alpha;
            beta_inv = 1d / beta;
        }

        public override ddouble PDF(ddouble x) {
            if (x < 0d || x > 1d) {
                return 0d;
            }

            if (Beta == 1d) {
                return Alpha * Pow(x, Alpha - 1d);
            }
            if (Alpha == 1d) {
                return Beta * Pow1p(-x, Beta - 1d);
            }

            if (x == 0d) {
                return Alpha < 1d ? PositiveInfinity : 0d;
            }

            if (x == 1d) {
                return Beta < 1d ? PositiveInfinity : 0d;
            }

            ddouble xa = Pow(x, Alpha);

            ddouble pdf = ab * xa * Pow1p(-xa, Beta - 1d) / x;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }
                if (x > 1d) {
                    return 1d;
                }

                ddouble cdf = 1d - Pow1p(-Pow(x, Alpha), Beta);

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }
                if (x > 1d) {
                    return 0d;
                }

                ddouble cdf = Pow1p(-Pow(x, Alpha), Beta);

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
            else {
                ddouble u = 1d - Pow(p, beta_inv);

                if (u <= 0d) {
                    return 0d;
                }
                if (u >= 1d) {
                    return 1d;
                }

                ddouble x = Pow(u, alpha_inv);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(1d - double.Pow(u, (double)beta_inv), (double)alpha_inv);

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean => Beta * Beta(1d + alpha_inv, Beta);

        public override ddouble Median => Pow1p(-Pow2(-beta_inv), alpha_inv);

        public override ddouble Mode =>
            Alpha >= 1d && Beta >= 1d && !(Alpha == 1d && Beta == 1d)
            ? Pow((Alpha - 1d) / (ab - 1d), alpha_inv)
            : NaN;

        public override ddouble Variance =>
            Beta * Beta(Beta, 1d + 2d * alpha_inv) - Square(Mean);

        public override ddouble Skewness {
            get {
                ddouble b1 = Beta(Beta, 1d + alpha_inv);
                ddouble b2 = Beta(Beta, 1d + 2d * alpha_inv);
                ddouble b3 = Beta(Beta, 1d + 3d * alpha_inv);

                return (2d * Cube(Beta * b1) - 3d * Beta * Beta * b1 * b2 + Beta * b3)
                     / ExMath.Pow3d2(Beta * b2 - Square(Beta * b1));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble b1 = Beta(Beta, 1d + alpha_inv);
                ddouble b2 = Beta(Beta, 1d + 2d * alpha_inv);
                ddouble b3 = Beta(Beta, 1d + 3d * alpha_inv);
                ddouble b4 = Beta(Beta, 1d + 4d * alpha_inv);

                return (-3d * Pow(Beta * b1, 4) + 6d * Cube(Beta) * b1 * b1 * b2 - 4d * Beta * Beta * b1 * b3 + Beta * b4)
                    / Square(Beta * b2 - Square(Beta * b1)) - 3d;
            }
        }

        public override ddouble Entropy =>
            (1d - beta_inv) + (1d - alpha_inv) * (Digamma(Beta + 1d) + EulerGamma) - Log(ab);

        public static (KumaraswamyDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (KumaraswamyDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble alpha = t.u / (1d - t.u);
                    ddouble beta = t.v / (1d - t.v);

                    try {
                        KumaraswamyDistribution dist = new(alpha, beta);
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
                KumaraswamyDistribution dist = new(alpha, beta);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(KumaraswamyDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := alpha * beta * x^(alpha - 1) * (1 - x^alpha)^(beta - 1)";
    }
}
