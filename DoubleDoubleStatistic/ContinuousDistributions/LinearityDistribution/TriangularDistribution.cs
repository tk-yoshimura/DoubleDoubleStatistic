using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
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
    public class TriangularDistribution : LinearityDistribution<TriangularDistribution>,
        IAdditionOperators<TriangularDistribution, ddouble, TriangularDistribution>,
        ISubtractionOperators<TriangularDistribution, ddouble, TriangularDistribution>,
        IMultiplyOperators<TriangularDistribution, ddouble, TriangularDistribution>,
        IDivisionOperators<TriangularDistribution, ddouble, TriangularDistribution>,
        IFittableContinuousDistribution<TriangularDistribution> {

        public ddouble A { get; }
        public ddouble B { get; }
        public ddouble C { get; }

        private readonly ddouble ab, ac, cb, abxac, abxcb, p_thr;

        public TriangularDistribution() : this(0d, 1d) { }

        public TriangularDistribution(ddouble a, ddouble b) : this(a, b, c: (a + b) * 0.5d) { }

        public TriangularDistribution(ddouble a, ddouble b, ddouble c) {
            ValidateLocation(a);
            ValidateLocation(b);
            ValidateLocation(c);

            if (!(a < c) || !(c < b) || !IsFinite(b - a)) {
                throw new ArgumentException($"Invalid location parameter. {nameof(a)} < {nameof(c)} < {nameof(b)}");
            }

            A = a;
            B = b;
            C = c;
            ab = b - a;
            ac = c - a;
            cb = b - c;
            abxac = ab * ac;
            abxcb = ab * cb;
            p_thr = ac / ab;
        }

        public override ddouble PDF(ddouble x) {
            if (x < A || x > B) {
                return 0d;
            }

            ddouble pdf = (x <= C)
                ? 2d * (x - A) / abxac
                : 2d * (B - x) / abxcb;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < A) {
                    return 0d;
                }

                if (x > B) {
                    return 1d;
                }

                ddouble cdf = (x <= C)
                    ? Square(x - A) / abxac
                    : 1d - Square(B - x) / abxcb;

                return cdf;
            }
            else {
                if (x < A) {
                    return 1d;
                }

                if (x > B) {
                    return 0d;
                }

                ddouble cdf = (x <= C)
                    ? 1d - Square(x - A) / abxac
                    : Square(B - x) / abxcb;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = (p < p_thr)
                    ? A + Sqrt(abxac * p)
                    : B - Sqrt(abxcb * (1d - p));

                return x;
            }
            else {
                ddouble x = ((1d - p) < p_thr)
                    ? A + Sqrt(abxac * (1d - p))
                    : B - Sqrt(abxcb * p);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniform();

            double v = (u < (double)p_thr)
                ? (double)A + double.Sqrt((double)abxac * u)
                : (double)B - double.Sqrt((double)abxcb * (1d - u));

            return v;
        }

        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mean => (A + B + C) / 3d;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => C;

        public override ddouble Variance => (A * A + B * B + C * C - A * B - A * C - B * C) / 18d;

        public override ddouble Skewness =>
            Sqrt2 * (2d * A - B - C) * (A - 2d * B + C) * (A + B - 2d * C)
            / (5d * ExMath.Pow3d2(A * A + B * B + C * C - A * B - A * C - B * C));

        public override ddouble Kurtosis => -(ddouble)3 / 5;

        public override ddouble Entropy => 0.5d + Log(ab * 0.5d);

        public static TriangularDistribution operator +(TriangularDistribution dist, ddouble s) {
            return new(dist.A + s, dist.B + s, dist.C + s);
        }

        public static TriangularDistribution operator -(TriangularDistribution dist, ddouble s) {
            return new(dist.A - s, dist.B - s, dist.C - s);
        }

        public static TriangularDistribution operator *(TriangularDistribution dist, ddouble k) {
            return new(dist.A * k, dist.B * k, dist.C * k);
        }

        public static TriangularDistribution operator /(TriangularDistribution dist, ddouble k) {
            return new(dist.A / k, dist.B / k, dist.C / k);
        }

        public static (TriangularDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (TriangularDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble c = GridMinimizeSearch1D.Search(
                c => {
                    try {
                        TriangularDistribution dist = new(0d, 1d, c);
                        return QuantileLinearFitter<TriangularDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (0.0001d, 0.9999d), iter: 32
            );

            try {
                TriangularDistribution dist = new(0d, 1d, c);

                return QuantileLinearFitter<TriangularDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(TriangularDistribution).Name}[a={A},b={B},c={C}]";
        }

        public override string Formula => "p(x; a, b, c) := if x < c then 2 * (x - a) / ((b - a) * (c - a)) else 2 * (b - x) / ((b - a) * (b - c))";
    }
}
