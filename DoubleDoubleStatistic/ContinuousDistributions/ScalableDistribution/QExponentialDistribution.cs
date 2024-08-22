using DoubleDouble;
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
    public class QExponentialDistribution : ScalableDistribution<QExponentialDistribution>,
        IMultiplyOperators<QExponentialDistribution, ddouble, QExponentialDistribution>,
        IDivisionOperators<QExponentialDistribution, ddouble, QExponentialDistribution>,
        IFittableContinuousDistribution<QExponentialDistribution> {

        public ddouble Q { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_power, cdf_power, quantile_power, pdf_norm, theta_inv, q12;

        public QExponentialDistribution(ddouble q) : this(q, theta: 1d) { }

        public QExponentialDistribution(ddouble q, ddouble theta) {
            ParamAssert.ValidateShape(nameof(q), q >= 0 && q < 2);
            ParamAssert.ValidateScale(nameof(theta), ParamAssert.IsFinitePositive(theta));

            Q = q;
            Theta = theta;

            theta_inv = 1d / theta;

            pdf_norm = (2d - Q) * theta_inv;
            pdf_power = 1d / (1d - q);
            cdf_power = (Q - 2d) / (Q - 1d);
            quantile_power = (Q - 1d) / (Q - 2d);
            q12 = (Q - 1d) * (Q - 2d);
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (IsNegative(u) || IsInfinity(u)) {
                return 0d;
            }

            if (IsFinite(pdf_power)) {
                ddouble pdf = pdf_norm * Pow(Max(0d, 1d + (Q - 1d) * u), pdf_power);

                if (IsNaN(pdf)) {
                    return 0d;
                }

                return pdf;
            }
            else {
                ddouble pdf = pdf_norm * Exp(-u);

                return pdf;
            }
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

                if (IsFinite(cdf_power)) {
                    ddouble cdf = 1d - Pow(Max(0d, 1d + (Q - 1d) * u), cdf_power);

                    cdf = Max(0d, cdf);

                    return cdf;
                }
                else {
                    ddouble cdf = -Expm1(-u);

                    return cdf;
                }
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                if (IsFinite(cdf_power)) {
                    ddouble cdf = Pow(Max(0d, 1d + (Q - 1d) * u), cdf_power);

                    cdf = Min(1d, cdf);

                    return cdf;
                }
                else {
                    ddouble cdf = Exp(-u);

                    return cdf;
                }
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble u;

            if (interval == Interval.Lower) {
                if (p >= 1d) {
                    return Support.max;
                }

                if (Q != 1d) {
                    u = (Pow1p(-p, quantile_power) - 1d) / (Q - 1d);
                }
                else {
                    u = -Log1p(-p);
                }
            }
            else {
                if (p <= 0d) {
                    return Support.max;
                }

                if (Q != 1d) {
                    u = (Pow(p, quantile_power) - 1d) / (Q - 1d);
                }
                else {
                    u = -Log(p);
                }
            }

            if (IsNegative(u)) {
                return 0d;
            }

            ddouble x = u * Theta;
            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double r = (Q != 1d)
                ? (double.Pow(u, (double)quantile_power) - 1d) / ((double)Q - 1d)
                : -double.Log(u);

            double w = r * (double)Theta;

            return w;
        }

        public override (ddouble min, ddouble max) Support => Q < 1d
            ? (0d, Theta / (1d - Q))
            : (0d, PositiveInfinity);

        public override ddouble Mean => (Q < 1.5d)
            ? Theta / (3d - 2d * Q)
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => (Q * 3d < 4d)
            ? Theta * Theta * (Q - 2d) / (Square(2d * Q - 3d) * (3d * Q - 4d))
            : NaN;

        public override ddouble Skewness => (Q < 1.25d)
            ? 2d / (5d - 4d * Q) * Sqrt((3d * Q - 4d) / (Q - 2d))
            : NaN;

        public override ddouble Kurtosis => (Q * 5d < 6d)
            ? 6d * (6d + Q * (-20d + Q * (17d + Q * -4d))) / ((Q - 2d) * (4d * Q - 5d) * (5d * Q - 6d))
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 32768);

        public static QExponentialDistribution operator *(QExponentialDistribution dist, ddouble k) {
            return new(dist.Q, dist.Theta * k);
        }

        public static QExponentialDistribution operator /(QExponentialDistribution dist, ddouble k) {
            return new(dist.Q, dist.Theta / k);
        }

        public static (QExponentialDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (QExponentialDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble q = BisectionMinimizeSearch1D.Search(
                q => {
                    try {
                        QExponentialDistribution dist = new(q, 1d);
                        return QuantileScaleFitter<QExponentialDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (0, 1.9995d), iter: 32
            );

            try {
                QExponentialDistribution dist = new(q, 1d);

                return QuantileScaleFitter<QExponentialDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(QExponentialDistribution).Name}[q={Q},theta={Theta}]";
        }

        public override string Formula => "p(x; q, theta) := (2d - Q) * e_q(-u) / theta, u = x / theta";
    }
}
