﻿using DoubleDouble;
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
    public class FiskDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<FiskDistribution> {

        public ddouble C { get; }

        private readonly ddouble c_inv;

        public FiskDistribution(ddouble c) {
            ParamAssert.ValidateShape(nameof(c), ParamAssert.IsFinitePositive(c));

            C = c;
            c_inv = 1d / c;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble xc = Pow(x, C);

            if (xc <= 0d) {
                return C < 1d ? PositiveInfinity : C == 1d ? 1d : 0d;
            }
            if (IsPositiveInfinity(xc)) {
                return 0d;
            }

            ddouble pdf = C * (xc / (x * Square(xc + 1d)));

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble xc = Pow(x, C), xcp1 = xc + 1d;

            if (interval == Interval.Lower) {
                if (IsNegative(x) || xc <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 1d;
                }

                ddouble cdf = xc / (1d + xc);

                return cdf;
            }
            else {
                if (IsNegative(x) || xc <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 0d;
                }

                ddouble cdf = 1d / (1d + xc);

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

                ddouble x = Pow(p / (1d - p), c_inv);

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Pow((1d - p) / p, c_inv);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(u / (1d - u), (double)c_inv);

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Beta(1d - c_inv, 1d + c_inv);

        public override ddouble Median => 1d;

        public override ddouble Mode =>
            Pow((C - 1d) / (C + 1d), c_inv);

        public override ddouble Variance {
            get {
                ddouble mu1 = Beta(1d - c_inv, 1d + c_inv);
                ddouble mu2 = Beta((C - 2d) * c_inv, (C + 2d) * c_inv);

                return mu2 - mu1 * mu1;
            }
        }
        public override ddouble Skewness {
            get {
                ddouble mu1 = Beta(1d - c_inv, 1d + c_inv);
                ddouble mu2 = Beta((C - 2d) * c_inv, (C + 2d) * c_inv);
                ddouble mu3 = Beta((C - 3d) * c_inv, (C + 3d) * c_inv);

                return (2d * Cube(mu1) - 3d * mu1 * mu2 + mu3) / ExMath.Pow3d2(mu2 - mu1 * mu1);
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mu1 = Beta(1d - c_inv, 1d + c_inv);
                ddouble mu2 = Beta((C - 2d) * c_inv, (C + 2d) * c_inv);
                ddouble mu3 = Beta((C - 3d) * c_inv, (C + 3d) * c_inv);
                ddouble mu4 = Beta((C - 4d) * c_inv, (C + 4d) * c_inv);

                return (-3d * Square(Square(mu1)) + 6d * mu1 * mu1 * mu2 - 4d * mu1 * mu3 + mu4) / Square(mu2 - mu1 * mu1);
            }
        }

        public override ddouble Entropy {
            get {
                ddouble f(ddouble x) {
                    ddouble pdf = PDF(x);

                    if (pdf == 0d) {
                        return 0d;
                    }

                    ddouble y = -pdf * Log(pdf);

                    return y;
                }

                (ddouble value, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                    f, 0, PositiveInfinity, 1e-28, discontinue_eval_points: 16384
                );

                return value;
            }
        }

        public static (FiskDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (FiskDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble c = t / (1d - t);

                    try {
                        FiskDistribution dist = new(c);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble c = t / (1d - t);
                FiskDistribution dist = new(c);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(FiskDistribution).Name}[c={C}]";
        }

        public override string Formula => "p(x; c) := c * x^(c - 1) / (1 + x^2)^2";
    }
}
