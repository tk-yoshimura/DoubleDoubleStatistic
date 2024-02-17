using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class TriangularDistribution : LinearityDistribution<TriangularDistribution>,
        IAdditionOperators<TriangularDistribution, ddouble, TriangularDistribution>,
        IMultiplyOperators<TriangularDistribution, ddouble, TriangularDistribution> {

        public ddouble A { get; }
        public ddouble B { get; }
        public ddouble C { get; }

        private readonly ddouble ab, ac, cb, abxac, abxcb, p_thr;

        public TriangularDistribution(ddouble a, ddouble b) : this(a, b, c: Ldexp(a + b, -1)) { }

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
                ? 2 * (x - A) / abxac
                : 2 * (B - x) / abxcb;

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

        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mean => (A + B + C) / 3d;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => C;

        public override ddouble Variance => (A * A + B * B + C * C - A * B - A * C - B * C) / 18d;

        public override ddouble Skewness =>
            Sqrt2 * (2 * A - B - C) * (A - 2 * B - C) * (A + B - 2 * C)
            / (5d * Cube(Sqrt(A * A + B * B + C * C - A * B - A * C - B * C)));

        public override ddouble Kurtosis => -(ddouble)3 / 5;

        public override ddouble Entropy => 0.5d + Log(ab / 2);

        public static TriangularDistribution operator +(TriangularDistribution dist, ddouble s) {
            return new(s + dist.A, s + dist.B, s + dist.C);
        }

        public static TriangularDistribution operator *(TriangularDistribution dist, ddouble k) {
            return new(k * dist.A, k * dist.B, k * dist.C);
        }

        public override string ToString() {
            return $"{typeof(TriangularDistribution).Name}[a={A},b={B},c={C}]";
        }
    }
}
