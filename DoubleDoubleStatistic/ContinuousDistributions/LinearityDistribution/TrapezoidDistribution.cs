using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class TrapezoidDistribution : LinearityDistribution<TrapezoidDistribution>,
        IAdditionOperators<TrapezoidDistribution, ddouble, TrapezoidDistribution>,
        ISubtractionOperators<TrapezoidDistribution, ddouble, TrapezoidDistribution>,
        IMultiplyOperators<TrapezoidDistribution, ddouble, TrapezoidDistribution>,
        IDivisionOperators<TrapezoidDistribution, ddouble, TrapezoidDistribution> {

        public ddouble A { get; }
        public ddouble B { get; }
        public ddouble C { get; }
        public ddouble D { get; }

        private readonly ddouble ab, cd, r, bp, cp;

        public TrapezoidDistribution(ddouble a, ddouble b, ddouble c, ddouble d) {
            ParamAssert.ValidateLocation(nameof(a), IsFinite(a));
            ParamAssert.ValidateLocation(nameof(b), IsFinite(b));
            ParamAssert.ValidateLocation(nameof(c), IsFinite(c));
            ParamAssert.ValidateLocation(nameof(d), IsFinite(d));
            ParamAssert.ValidateLocation($"{nameof(a)},{nameof(b)},{nameof(c)},{nameof(d)}", (a < b) && (b < c) && (c < d) && IsFinite(d - a));

            A = a;
            B = b;
            C = c;
            D = d;
            ab = b - a;
            cd = d - c;
            r = 1d / (d + c - a - b);

            bp = CDF(B);
            cp = CDF(C);
        }

        public override ddouble PDF(ddouble x) {
            if (x < A || x > D) {
                return 0d;
            }

            ddouble pdf =
                (x <= B) ? 2d * r * (x - A) / ab :
                (x <= C) ? 2d * r :
                2d * r * (D - x) / cd;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < A) {
                    return 0d;
                }
                if (x > D) {
                    return 1d;
                }

                ddouble cdf =
                    (x <= B) ? r * Square(x - A) / ab :
                    (x <= C) ? r * (2d * x - A - B) :
                    1d - r * Square(D - x) / cd;

                return cdf;
            }
            else {
                if (x < A) {
                    return 1d;
                }
                if (x > D) {
                    return 0d;
                }

                ddouble cdf =
                    (x <= B) ? 1d - r * Square(x - A) / ab :
                    (x <= C) ? 1d - r * (2d * x - A - B) :
                    r * Square(D - x) / cd;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p, Interval.Lower);
            }

            if (p < bp) {
                ddouble x = A + Sqrt(ab * p / r);

                return x;
            }
            else if (p < cp) {
                ddouble x = ((A + B) * r + p) / (2d * r);

                return x;
            }
            else {
                p = 1d - p;
                ddouble x = D - Sqrt(cd * p / r);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniform();

            if (u < bp) {
                double x = (double)A + double.Sqrt((double)ab * u / (double)r);

                return x;
            }
            else if (u < cp) {
                double x = (((double)A + (double)B) * (double)r + u) / (2d * (double)r);

                return x;
            }
            else {
                u = 1d - u;
                double x = (double)D - double.Sqrt((double)cd * u / (double)r);

                return x;
            }
        }

        public override bool Symmetric => ab == cd;

        public override (ddouble min, ddouble max) Support => (A, D);

        public override ddouble Mean =>
            ((Cube(D) - Cube(C)) / cd - (Cube(B) - Cube(A)) / ab) / 3d * r;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => NaN;

        public override ddouble Variance =>
            ((Pow(D, 4) - Pow(C, 4)) / cd - (Pow(B, 4) - Pow(A, 4)) / ab) / 6d * r - Square(Mean);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 8192);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 8192);

        public override ddouble Entropy =>
            (D - C + B - A) / 2d * r - Log(r * 2d);

        public static TrapezoidDistribution operator +(TrapezoidDistribution dist, ddouble s) {
            return new(dist.A + s, dist.B + s, dist.C + s, dist.D + s);
        }

        public static TrapezoidDistribution operator -(TrapezoidDistribution dist, ddouble s) {
            return new(dist.A - s, dist.B - s, dist.C - s, dist.D - s);
        }

        public static TrapezoidDistribution operator *(TrapezoidDistribution dist, ddouble k) {
            return new(dist.A * k, dist.B * k, dist.C * k, dist.D * k);
        }

        public static TrapezoidDistribution operator /(TrapezoidDistribution dist, ddouble k) {
            return new(dist.A / k, dist.B / k, dist.C / k, dist.D / k);
        }

        public override string ToString() {
            return $"{typeof(TrapezoidDistribution).Name}[a={A},b={B},c={C},d={D}]";
        }

        public override string Formula => "p(x; a, b, c, d) := if x < c then 2 * (x - a) / ((b - a) * (c - a)) else 2 * (b - x) / ((b - a) * (b - c))";
    }
}
