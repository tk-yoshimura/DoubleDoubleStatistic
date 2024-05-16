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
    public class BurrDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<BurrDistribution> {

        public ddouble C { get; }
        public ddouble K { get; }

        private readonly ddouble ck, c_inv, k_inv;

        public BurrDistribution(ddouble c, ddouble k) {
            ParamAssert.ValidateShape(nameof(c), ParamAssert.IsFinitePositive(c));
            ParamAssert.ValidateShape(nameof(k), ParamAssert.IsFinitePositive(k));

            C = c;
            K = k;

            ck = c * k;
            c_inv = 1d / c;
            k_inv = 1d / k;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (C == 1d) {
                ddouble pdf = K / Pow(x + 1d, K + 1d);
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }
            else {
                ddouble xc = Pow(x, C);

                if (xc <= 0d) {
                    return C < 1d ? PositiveInfinity : 0d;
                }

                if (IsPositiveInfinity(xc)) {
                    return 0d;
                }

                ddouble pdf = ck * xc / (x * Pow(xc + 1d, K + 1d));
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble xc = Pow(x, C), xcp1 = 1d + xc;

            if (interval == Interval.Lower) {
                if (IsNegative(x) || xc <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 1d;
                }

                ddouble cdf = Max(0d, 1d - Pow(xcp1, -K));

                return cdf;
            }
            else {
                if (IsNegative(x) || xc <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 0d;
                }

                ddouble cdf = Min(1d, Pow(xcp1, -K));

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

                ddouble x = Pow(Pow(1d / (1d - p), k_inv) - 1d, c_inv);

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Pow(Pow(1d / p, k_inv) - 1d, c_inv);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(double.Pow(1d / u, (double)k_inv) - 1d, (double)c_inv);

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            K * Beta(K - c_inv, 1d + c_inv);

        public override ddouble Median =>
            Pow(Pow2(k_inv) - 1d, c_inv);

        public override ddouble Mode =>
            Pow((C - 1d) / (K * C + 1d), c_inv);

        public override ddouble Variance {
            get {
                ddouble mu1 = K * Beta(K - c_inv, 1d + c_inv);
                ddouble mu2 = K * Beta((C * K - 2d) * c_inv, (C + 2d) * c_inv);

                return mu2 - mu1 * mu1;
            }
        }
        public override ddouble Skewness {
            get {
                ddouble mu1 = K * Beta(K - c_inv, 1d + c_inv);
                ddouble mu2 = K * Beta((C * K - 2d) * c_inv, (C + 2d) * c_inv);
                ddouble mu3 = K * Beta((C * K - 3d) * c_inv, (C + 3d) * c_inv);

                return (2d * Cube(mu1) - 3d * mu1 * mu2 + mu3) / ExMath.Pow3d2(mu2 - mu1 * mu1);
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mu1 = K * Beta(K - c_inv, 1d + c_inv);
                ddouble mu2 = K * Beta((C * K - 2d) * c_inv, (C + 2d) * c_inv);
                ddouble mu3 = K * Beta((C * K - 3d) * c_inv, (C + 3d) * c_inv);
                ddouble mu4 = K * Beta((C * K - 4d) * c_inv, (C + 4d) * c_inv);

                return (-3d * Square(Square(mu1)) + 6d * mu1 * mu1 * mu2 - 4d * mu1 * mu3 + mu4) / Square(mu2 - mu1 * mu1) - 3d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (BurrDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BurrDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = GridMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble c = t.u / (1d - t.u);
                    ddouble k = t.v / (1d - t.v);

                    try {
                        BurrDistribution dist = new(c, k);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-10d, 1e-10d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble c = u / (1d - u);
                ddouble k = v / (1d - v);
                BurrDistribution dist = new(c, k);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(BurrDistribution).Name}[c={C},k={K}]";
        }

        public override string Formula => "p(x; c, k) := c * k * x^(c - 1) / (1 + x^c)^(k + 1)";
    }
}
