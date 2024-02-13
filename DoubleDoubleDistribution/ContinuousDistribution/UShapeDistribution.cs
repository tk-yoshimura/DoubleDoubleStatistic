using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class UShapeDistribution : ContinuousDistribution,
        IAdditionOperators<UShapeDistribution, ddouble, UShapeDistribution>,
        IMultiplyOperators<UShapeDistribution, ddouble, UShapeDistribution> {

        public ddouble A { get; }
        public ddouble B { get; }

        private readonly ddouble alpha, beta, c, range;

        public UShapeDistribution(ddouble a, ddouble b) {
            ValidateLocation(a);
            ValidateLocation(b);

            if (!(a < b) || !IsFinite(b - a)) {
                throw new ArgumentException($"Invalid location parameter. {nameof(a)} < {nameof(b)}");
            }

            A = a;
            B = b;
            range = b - a;

            alpha = 12d / Cube(range);
            beta = Ldexp(A + B, -1);
            c = Cube(beta - A);
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
                ddouble x = beta + Cbrt(3d * p - alpha * c) / Cbrt(alpha);

                x = Clamp(x, A, B);

                return x;
            }
            else {
                ddouble x = beta - Cbrt(3d * p - alpha * c) / Cbrt(alpha);

                x = Clamp(x, A, B);

                return x;
            }
        }

        public override bool Scalable => true;
        public override bool Shiftable => true;
        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mean => Ldexp(A + B, -1);

        public override ddouble Median => Ldexp(A + B, -1);

        public override ddouble Mode => NaN;

        public override ddouble Variance =>
            Square(range) * 3d / 20d;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis =>
            Square(Square(range)) * 3d / 112d;

        public override ddouble Entropy {
            get {
                ddouble r = 3 * Log(3 / range);

                return -(Cube(Sqrt(range)) * Exp(r / 2) * (r - 2)) / (3 * Cube(Sqrt(3)));
            }
        }

        public static UShapeDistribution operator +(UShapeDistribution dist, ddouble s) {
            return new(s + dist.A, s + dist.B);
        }

        public static UShapeDistribution operator *(UShapeDistribution dist, ddouble k) {
            return new(k * dist.A, k * dist.B);
        }

        public override string ToString() {
            return $"{typeof(UShapeDistribution).Name}[a={A},b={B}]";
        }
    }
}
