using DoubleDouble;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class UniformDistribution : LinearityDistribution<UniformDistribution>,
        IAdditionOperators<UniformDistribution, ddouble, UniformDistribution>,
        ISubtractionOperators<UniformDistribution, ddouble, UniformDistribution>,
        IMultiplyOperators<UniformDistribution, ddouble, UniformDistribution>,
        IDivisionOperators<UniformDistribution, ddouble, UniformDistribution> {

        public ddouble A { get; }
        public ddouble B { get; }

        private readonly ddouble pdf_norm, range;

        public UniformDistribution() : this(0d, 1d) { }

        public UniformDistribution(ddouble a, ddouble b) {
            ValidateLocation(a);
            ValidateLocation(b);

            if (!(a < b) || !IsFinite(b - a)) {
                throw new ArgumentException($"Invalid location parameter. {nameof(a)} < {nameof(b)}");
            }

            A = a;
            B = b;

            range = b - a;
            pdf_norm = 1d / range;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x < A || x > B) {
                return 0d;
            }

            ddouble pdf = pdf_norm;

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

                ddouble cdf = (x - A) * pdf_norm;

                return cdf;
            }
            else {
                if (x < A) {
                    return 1d;
                }

                if (x > B) {
                    return 0d;
                }

                ddouble cdf = (B - x) * pdf_norm;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = A + p * range;

                return x;
            }
            else {
                ddouble x = B - p * range;

                return x;
            }
        }

        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mode => NaN;

        public override ddouble Mean => (A + B) * 0.5d;

        public override ddouble Median => (A + B) * 0.5d;

        public override ddouble Variance => Square(range) / 12d;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => -(ddouble)6 / 5;

        public override ddouble Entropy => Log(range);

        public static UniformDistribution operator +(UniformDistribution dist, ddouble s) {
            return new(dist.A + s, dist.B + s);
        }

        public static UniformDistribution operator -(UniformDistribution dist, ddouble s) {
            return new(dist.A - s, dist.B - s);
        }

        public static UniformDistribution operator *(UniformDistribution dist, ddouble k) {
            return new(dist.A * k, dist.B * k);
        }

        public static UniformDistribution operator /(UniformDistribution dist, ddouble k) {
            return new(dist.A / k, dist.B / k);
        }

        public override string ToString() {
            return $"{typeof(UniformDistribution).Name}[a={A},b={B}]";
        }

        public override string Formula => "p(x; a, b) := 1 / (b - a)";
    }
}
