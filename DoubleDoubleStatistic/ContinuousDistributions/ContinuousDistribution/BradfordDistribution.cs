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
    public class BradfordDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<BradfordDistribution> {

        public ddouble C { get; }

        private readonly ddouble pdf_norm, c_inv;

        public BradfordDistribution(ddouble c) {
            ParamAssert.ValidateShape(nameof(c), ParamAssert.IsFinitePositive(c));

            C = c;
            c_inv = 1d / c;
            pdf_norm = c / Log1p(c);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || x > 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm / (1d + C * x);

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

                ddouble cdf = Log1p(C * x) * pdf_norm * c_inv;

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                if (x > 1d) {
                    return 0d;
                }

                ddouble cdf = 1d - Log1p(C * x) * pdf_norm * c_inv;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p);
            }

            ddouble x = Expm1(p * C / pdf_norm) * c_inv;

            x = Clamp(x, 0d, 1d);

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = (double.Exp(u * (double)C / (double)pdf_norm) - 1d) * (double)c_inv;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean {
            get {
                ddouble k = Log1p(C);
                return (C - k) / (C * k);
            }
        }

        public override ddouble Median => Expm1(C / (pdf_norm * 2d)) * c_inv;

        public override ddouble Mode => 0d;

        public override ddouble Variance {
            get {
                ddouble k = Log1p(C);
                return ((C + 2d) * k - 2d * C) / (2d * C * k * k);
            }
        }

        public override ddouble Skewness {
            get {
                ddouble k = Log1p(C);
                return Sqrt2 * (2d * k * k * (3d + C * (3d + C)) - 9d * C * k * (C + 2d) + 12d * C * C) /
                    (Sqrt(C * (2d * k + C * (k - 2d))) * (6d * k + 3d * C * (k - 2d)));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble k = Log1p(C);
                return (12d * k * k * k
                    + 6d * C * k * k * (3d * k - 14d)
                    + 12d * C * C * k * (k - 4d) * (k - 3d)
                    + C * C * C * (k - 3d) * (3d * k * k - 16d * k + 24d))
                    / (3d * C * Square(2d * k + C * (k - 2d)));
            }
        }

        public override ddouble Entropy {
            get {
                ddouble k = Log1p(C);
                return k * 0.5d - Log(C / k);
            }
        }

        public static (BradfordDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BradfordDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble c = t / (1d - t);

                    try {
                        BradfordDistribution dist = new(c);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble c = t / (1d - t);
                BradfordDistribution dist = new(c);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(BradfordDistribution).Name}[c={C}]";
        }

        public override string Formula => "p(x; c) := c / (log(1 + c) * (1 + c * x))";
    }
}
