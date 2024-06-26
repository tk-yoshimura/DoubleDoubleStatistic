﻿using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class SnedecorFDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<SnedecorFDistribution> {

        public ddouble N { get; }
        public ddouble M { get; }

        private readonly ddouble pdf_lognorm;

        private readonly ChiSquareDistribution randam_gen_chisq_dist_n, randam_gen_chisq_dist_m;

        public SnedecorFDistribution(ddouble n, ddouble m) {
            ParamAssert.ValidateShape(nameof(n), ParamAssert.IsFinitePositive(n));
            ParamAssert.ValidateShape(nameof(m), ParamAssert.IsFinitePositive(m));

            N = n;
            M = m;

            pdf_lognorm = m * Log2(m) * 0.5d - LogBeta(n * 0.5d, m * 0.5d) * LbE;

            randam_gen_chisq_dist_n = new(n);
            randam_gen_chisq_dist_m = new(m);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x)) {
                return 0d;
            }

            ddouble u = N * x;

            if (u <= 0d) {
                return N < 2d ? PositiveInfinity : N == 2d ? 1d : 0d;
            }

            if (IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble pdf = Pow2((N * Log2(u) - (N + M) * Log2(u + M)) * 0.5d + pdf_lognorm) / x;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = N * x, upm = u + M;

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(upm)) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(u / upm, N * 0.5d, M * 0.5d);

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(upm)) {
                    return 0d;
                }

                ddouble cdf = IncompleteBetaRegularized(M / upm, M * 0.5d, N * 0.5d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = InverseIncompleteBeta(p, N * 0.5d, M * 0.5d);
                ddouble x = M * u / (N * (1d - u));

                return x;
            }
            else {
                ddouble u = InverseIncompleteBeta(p, M * 0.5d, N * 0.5d);
                ddouble x = M * (1d - u) / (N * u);

                return x;
            }
        }

        public override double Sample(Random random) {
            double c1 = randam_gen_chisq_dist_n.Sample(random);
            double c2 = randam_gen_chisq_dist_m.Sample(random);

            return (c1 * (double)M) / double.Max(c2 * (double)N, double.Epsilon);
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => M > 2d
            ? M / (M - 2d)
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => N > 2d
            ? (N - 2d) * M / (N * (M + 2d))
            : NaN;

        public override ddouble Variance => M > 4d
            ? 2d * M * M * (N + M - 2d) / (N * (M - 2d) * (M - 2d) * (M - 4d))
            : NaN;

        public override ddouble Skewness => M > 6d
            ? (2d * N + M - 2d) * Sqrt(8d * (M - 4d)) / ((M - 6d) * Sqrt(N * (N + M - 2d)))
            : NaN;

        public override ddouble Kurtosis => M > 8d
            ? 12d * (N * (5d * M - 22d) * (N + M - 2d) + (M - 4d) * (M - 2d) * (M - 2d))
                / (N * (M - 6d) * (M - 8d) * (N + M - 2d))
            : NaN;

        public override ddouble Entropy =>
            LogGamma(N * 0.5d) + LogGamma(M * 0.5d) - LogGamma((N + M) * 0.5d)
            + (1d - N * 0.5d) * Digamma(N * 0.5d) - (1d + M * 0.5d) * Digamma(M * 0.5d)
            + (N + M) * 0.5d * Digamma((N + M) * 0.5d) - Log(N / M);

        public static (SnedecorFDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (SnedecorFDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble n = t.u / (1d - t.u);
                    ddouble m = t.v / (1d - t.v);

                    try {
                        SnedecorFDistribution dist = new(n, m);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-4d, 1e-4d), (1000d / 1001d, 1000d / 1001d)), iter: 64
            );

            try {
                ddouble n = u / (1d - u);
                ddouble m = v / (1d - v);
                SnedecorFDistribution dist = new(n, m);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(SnedecorFDistribution).Name}[n={N},m={M}]";
        }

        public override string Formula => "p(x; n, m) := x^(n / 2 - 1) * (m + n * x)^((- m - n) / 2) * m^(m / 2) * n^(n / 2) / beta(n / 2, m / 2)";
    }
}
