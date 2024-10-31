using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class ArcsineDistribution : ContinuousDistribution {

        public ArcsineDistribution() { }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }
            if (IsNegative(x) || x > 1d) {
                return 0d;
            }

            ddouble pdf = 1d / (Pi * Sqrt(x * (1d - x)));

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }

                if (x > 1d) {
                    return 1d;
                }

                ddouble cdf = 2d * RcpPi * Asin(Sqrt(x));

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                if (x > 1d) {
                    return 0d;
                }

                ddouble cdf = 2d * RcpPi * Asin(Sqrt(1d - x));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p);
            }

            ddouble x = Square(SinPi(p * 0.5d));

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.SinPi(u * 0.5d);
            double w = v * v;

            return w;
        }

        public override bool Symmetric => true;

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean => 0.5d;

        public override ddouble Median => 0.5d;

        public override ddouble Mode => NaN;

        public override ddouble Variance => 0.125d;

        public override ddouble Skewness => 0;

        public override ddouble Kurtosis => -1.5d;

        public override ddouble Entropy => Log(Pi * 0.25d);

        public override string ToString() {
            return $"{typeof(ArcsineDistribution).Name}[]";
        }

        public override string Formula => "p(x) := 1 / (pi * sqrt(x * (1 - x)))";
    }
}
