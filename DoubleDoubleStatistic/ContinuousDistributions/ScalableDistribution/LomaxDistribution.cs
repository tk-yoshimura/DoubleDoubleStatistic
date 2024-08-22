﻿using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LomaxDistribution : ScalableDistribution<LomaxDistribution>,
        IMultiplyOperators<LomaxDistribution, ddouble, LomaxDistribution>,
        IDivisionOperators<LomaxDistribution, ddouble, LomaxDistribution>,
        IFittableContinuousDistribution<LomaxDistribution> {

        public ddouble Alpha { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_norm, alpha_inv, theta_inv;

        public LomaxDistribution(ddouble alpha) : this(alpha, theta: 1d) { }

        public LomaxDistribution(ddouble alpha, ddouble theta) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateScale(nameof(theta), ParamAssert.IsFinitePositive(theta));

            Alpha = alpha;
            Theta = theta;

            pdf_norm = alpha / theta;
            alpha_inv = 1d / alpha;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble pdf = pdf_norm / Pow1p(u, Alpha + 1d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }
                if (IsPositiveInfinity(u)) {
                    return 1d;
                }

                ddouble cdf = 1d - 1d / Pow1p(u, Alpha);

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = 1d / Pow1p(u, Alpha);

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

            if (p == 0d) {
                return PositiveInfinity;
            }
            if (p == 1d) {
                return 0d;
            }

            ddouble x = Theta * (Pow(p, -alpha_inv) - 1d);

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = (double)Theta * (double.Pow(u, -(double)alpha_inv) - 1d);

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d)
            ? Theta / (Alpha - 1d)
            : PositiveInfinity;

        public override ddouble Median =>
            Theta * (Pow(0.5d, -alpha_inv) - 1d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => (Alpha > 2d)
            ? Alpha * Theta * Theta / (Square(Alpha - 1d) * (Alpha - 2d))
            : (Alpha > 1d) ? PositiveInfinity : NaN;

        public override ddouble Skewness => (Alpha > 3d)
            ? 2d * (Alpha + 1d) / (Alpha - 3d) * Sqrt((Alpha - 2d) / Alpha)
            : NaN;

        public override ddouble Kurtosis => (Alpha > 4d)
            ? 6d * (-2d + Alpha * (-6d + Alpha * (1d + Alpha))) / (Alpha * (Alpha - 3d) * (Alpha - 4d))
            : NaN;

        public override ddouble Entropy => 1d + Log(Theta / Alpha) + 1d / Alpha;

        public static LomaxDistribution operator *(LomaxDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Theta * k);
        }

        public static LomaxDistribution operator /(LomaxDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Theta / k);
        }

        public static (LomaxDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (LomaxDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble alpha = t / (1d - t);

                    try {
                        LomaxDistribution dist = new(alpha, 1d);
                        return QuantileScaleFitter<LomaxDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble alpha = t / (1d - t);
                LomaxDistribution dist = new(alpha, 1d);

                return QuantileScaleFitter<LomaxDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(LomaxDistribution).Name}[alpha={Alpha},theta={Theta}]";
        }

        public override string Formula => "p(x; alpha, theta) := (1 + u)^(-alpha - 1) * (alpha / theta), u = x / theta";
    }
}
