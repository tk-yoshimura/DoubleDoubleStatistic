using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class UQuadraticDistribution : LinearityDistribution<UQuadraticDistribution>,
        IAdditionOperators<UQuadraticDistribution, ddouble, UQuadraticDistribution>,
        ISubtractionOperators<UQuadraticDistribution, ddouble, UQuadraticDistribution>,
        IMultiplyOperators<UQuadraticDistribution, ddouble, UQuadraticDistribution>,
        IDivisionOperators<UQuadraticDistribution, ddouble, UQuadraticDistribution>,
        IFittableDistribution<UQuadraticDistribution> {

        public ddouble A { get; }
        public ddouble B { get; }

        private readonly ddouble alpha, beta, c, range;

        public UQuadraticDistribution() : this(0d, 1d) { }

        public UQuadraticDistribution(ddouble a, ddouble b) {
            ValidateLocation(a);
            ValidateLocation(b);

            if (!(a < b) || !IsFinite(b - a)) {
                throw new ArgumentException($"Invalid location parameter. {nameof(a)} < {nameof(b)}");
            }

            A = a;
            B = b;
            range = b - a;

            alpha = 12d / Cube(range);
            beta = (a + b) * 0.5d;
            c = Cube(beta - a);
        }

        public override ddouble PDF(ddouble x) {
            if (x < A || x > B) {
                return 0d;
            }

            ddouble pdf = alpha * Square(x - beta);

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

                ddouble cdf = alpha * (Cube(x - beta) + c) / 3d;

                return cdf;
            }
            else {
                if (x < A) {
                    return 1d;
                }

                if (x > B) {
                    return 0d;
                }

                ddouble cdf = alpha * (Cube(beta - x) + c) / 3d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = beta + Cbrt(3d * p / alpha - c);

                x = Clamp(x, A, B);

                return x;
            }
            else {
                ddouble x = beta - Cbrt(3d * p / alpha - c);

                x = Clamp(x, A, B);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniform();

            double v = (double)beta + double.Cbrt(3d * u / (double)alpha - (double)c);

            return v;
        }

        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mean => (A + B) * 0.5d;

        public override ddouble Median => (A + B) * 0.5d;

        public override ddouble Mode => NaN;

        public override ddouble Variance =>
            Square(range) * 3d / 20d;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis =>
            -(ddouble)38d / 21d;

        public override ddouble Entropy {
            get {
                ddouble r = 3d * Log(3d / range);

                return -(ExMath.Pow3d2(range) * Exp(r * 0.5d) * (r - 2d)) / (3d * ExMath.Pow3d2(3d));
            }
        }

        public static UQuadraticDistribution operator +(UQuadraticDistribution dist, ddouble s) {
            return new(dist.A + s, dist.B + s);
        }

        public static UQuadraticDistribution operator -(UQuadraticDistribution dist, ddouble s) {
            return new(dist.A - s, dist.B - s);
        }

        public static UQuadraticDistribution operator *(UQuadraticDistribution dist, ddouble k) {
            return new(dist.A * k, dist.B * k);
        }

        public static UQuadraticDistribution operator /(UQuadraticDistribution dist, ddouble k) {
            return new(dist.A / k, dist.B / k);
        }

        public static (UQuadraticDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (UQuadraticDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileLinearFitter<UQuadraticDistribution>.Fit(new UQuadraticDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(UQuadraticDistribution).Name}[a={A},b={B}]";
        }

        public override string Formula => "p(x; a, b) := 12 / (b - a)^3 * (x - (a + b) / 2)^2";
    }
}
