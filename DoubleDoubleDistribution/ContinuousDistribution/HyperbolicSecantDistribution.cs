using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class HyperbolicSecantDistribution : ContinuousDistribution {

        public HyperbolicSecantDistribution() { }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = 0.5d / Cosh(x * PI * 0.5d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble cdf = 2 * RcpPI * Atan(Exp(-Abs(x) * PI * 0.5d));
            cdf = Max(cdf, 0d);

            cdf = (interval != Interval.Lower ^ IsNegative(x)) ? cdf : 1d - cdf;

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble x = Log(TanPI(p * 0.5d)) * RcpPI * 2d;

            x = (interval == Interval.Lower) ? x : -x;

            return x;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => 0d;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Variance => 1d;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => 2d;

        public override ddouble Entropy => Log(4d);

        public override string ToString() {
            return $"{typeof(HyperbolicSecantDistribution).Name}[]";
        }
    }
}
