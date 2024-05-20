using DoubleDouble;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BatesDistribution : ContinuousDistribution {

        public int N => dist.N;

        private readonly IrwinHallDistribution dist;

        public BatesDistribution(int n) {
            dist = new IrwinHallDistribution(n);
        }

        public override ddouble PDF(ddouble x) {
            return dist.PDF(x * N) * N;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            return dist.CDF(x * N, interval);
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            return dist.Quantile(p, interval) / N;
        }

        public override double Sample(Random random) {
            double u = dist.Sample(random) / N;

            return u;
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean => 0.5d;

        public override ddouble Median => 0.5d;

        public override ddouble Mode => (N > 1) ? 0.5d : NaN;

        public override ddouble Variance => 1d / (ddouble)(12 * N);

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => -(ddouble)6 / (5 * N);

        public override ddouble Entropy => dist.Entropy - Log(N);

        public override string ToString() {
            return $"{typeof(BatesDistribution).Name}[n={N}]";
        }

        public override string Formula => "p(x; n) := irwinhall(x * N; N) * N";
    }
}
