﻿using DoubleDouble;
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
    public class RiceDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<RiceDistribution> {

        public ddouble Nu { get; }

        private const int cache_samples = 512;
        private CDFSegmentCache? cdf_cache;
        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null;

        public RiceDistribution(ddouble nu) {
            ParamAssert.ValidateShape(nameof(nu), nu >= 0d && IsFinite(nu));

            Nu = nu;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble pdf = x * Exp(-Square(x - Nu) * 0.5d) * BesselI(0, x * Nu, scale: true);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d) {
                return interval == Interval.Lower ? 0d : 1d;
            }

            ddouble mode = Mode;

            if (IsPositiveInfinity(x + mode)) {
                return interval == Interval.Lower ? 1d : 0d;
            }

            this.cdf_cache ??= new CDFSegmentCache(0d, 1d,
                t => t > 0d ? PDF((1d - t) / t * mode) / (t * t) : 0d,
                samples: cache_samples
            );

            ddouble t = mode / (x + mode);

            ddouble cdf = (interval == Interval.Lower)
                ? Min(1d, cdf_cache.Upper(t) * mode)
                : Min(1d, cdf_cache.Lower(t) * mode);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            ddouble mode = Mode;

            this.cdf_cache ??= new CDFSegmentCache(0d, 1d,
                t => t > 0d ? PDF((1d - t) / t * mode) / (t * t) : 0d,
                samples: cache_samples
            );

            if (interval == Interval.Lower) {
                this.quantile_lower_builder ??= new QuantileBuilder(
                    1d, 0d,
                    new ReadOnlyCollection<ddouble>(cdf_cache.UpperSegmentTable.Reverse().ToArray()),
                    cdf_cache.Samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_lower_builder.Estimate(p / mode);

                ddouble x = (1d - t) / t * mode;
                ddouble x0 = (1d - t0) / t0 * mode;
                ddouble x1 = (1d - t1) / t1 * mode;

                for (int i = 0; i < 8; i++) {
                    ddouble y = CDF(x, Interval.Lower), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x - dx, x0, x1);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                this.quantile_upper_builder ??= new QuantileBuilder(
                    0d, 1d,
                    cdf_cache.LowerSegmentTable,
                    cdf_cache.Samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p / mode);

                ddouble x = (1d - t) / t * mode;
                ddouble x0 = (1d - t0) / t0 * mode;
                ddouble x1 = (1d - t1) / t1 * mode;

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

                return x;
            }
        }

        public override double Sample(Random random) {
            double theta = 2d * random.NextUniformOpenInterval1();
            (double s, double c) = double.SinCosPi(theta);

            (double u0, double u1) = random.NextGaussianX2();

            double x = c * (double)Nu + u0;
            double y = s * (double)Nu + u1;

            double r = double.Hypot(x, y);

            return r;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble x = Nu + 1d / (Nu + 1d);

                for (int i = 0; i < 32; i++) {
                    ddouble b0 = BesselI(0, x * Nu, scale: true), b1 = BesselI(1, x * Nu, scale: true);

                    ddouble dx = (x * Nu * b1 - (x * x - 1d) * b0) / (x * (x * x + Nu * Nu - 3d) * b0 + Nu * (1d - 2d * x * x) * b1);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x -= dx;

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                mode ??= x;

                return x;
            }
        }

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static (RiceDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (RiceDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble nu = t / (1d - t);

                    try {
                        RiceDistribution dist = new(nu);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (0d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble nu = t / (1d - t);
                RiceDistribution dist = new(nu);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(RiceDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := x * exp(-(x^2 + v^2) / 2) * bessel_i(0, x * nu)";
    }
}
