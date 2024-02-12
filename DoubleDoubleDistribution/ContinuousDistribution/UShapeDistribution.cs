using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class UShapeDistribution : ContinuousDistribution {

        public ddouble Min { get; }
        public ddouble Max { get; }
        public ddouble Range { get; }

        private readonly ddouble alpha, beta, c;

        public UShapeDistribution(ddouble min, ddouble max) {
            ValidateLocation(min);
            ValidateLocation(max);

            if (!(min < max)) {
                throw new ArgumentException($"Invalid location parameter. {min} < {max}");
            }

            Min = min;
            Max = max;
            Range = max - min;

            alpha = 12d / Cube(Range);
            beta = Ldexp(Min + Max, -1);
            c = Cube(beta - Min);
        }

        public override ddouble PDF(ddouble x) {
            if (x < Min || x > Max) {
                return 0d;
            }

            ddouble pdf = alpha * Square(x - beta);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < Min) {
                    return 0d;
                }

                if (x > Max) {
                    return 1d;
                }

                ddouble cdf = alpha * (Cube(x - beta) + c) / 3d;

                return cdf;
            }
            else {
                if (x < Min) {
                    return 1d;
                }

                if (x > Max) {
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
                ddouble quantile = beta + Cbrt(3d * p - alpha * c) / Cbrt(alpha);

                quantile = Clamp(quantile, Min, Max);

                return quantile;
            }
            else {
                ddouble quantile = beta - Cbrt(3d * p - alpha * c) / Cbrt(alpha);

                quantile = Clamp(quantile, Min, Max);

                return quantile;
            }
        }

        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (Min, Max);

        public override ddouble Mean => Ldexp(Min + Max, -1);
        public override ddouble Median => Ldexp(Min + Max, -1);
        public override ddouble Mode => NaN;

        public override ddouble Variance => Square(Range) * 3d / 20d;
        public override ddouble Skewness => 0d;
        public override ddouble Kurtosis => Square(Square(Range)) * 3d / 112d;

        public override ddouble Entropy {
            get {
                ddouble r = 3 * Log(3 / Range);

                return -(Cube(Sqrt(Range)) * Exp(r / 2) * (r - 2)) / (3 * Cube(Sqrt(3)));
            }
        }

        public override string ToString() {
            return $"{typeof(UShapeDistribution).Name}[min={Min},max={Max}]";
        }
    }
}
