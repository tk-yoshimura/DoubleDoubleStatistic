using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class PascalDistribution : DiscreteDistribution,
        IFittableDiscreteDistribution<PascalDistribution> {

        public int N => dist.N;
        public ddouble P => dist.P;

        private readonly NegativeBinomialDistribution dist;

        public PascalDistribution(int n, ddouble p) {
            ParamAssert.ValidateShape(nameof(n), n > 0);
            ParamAssert.ValidateShape(nameof(p), p >= 0d && p <= 1d);

            dist = new(n, p);
        }

        public override ddouble PMF(int k) {
            return dist.PMF(k - N);
        }

        public override int Sample(Random random) {
            return checked(dist.Sample(random) + N);
        }

        public override ddouble Mean => N / P;

        public override ddouble Variance => dist.Variance;

        public override ddouble Skewness => dist.Skewness;

        public override ddouble Kurtosis => dist.Kurtosis;

        public override ddouble Entropy => dist.Entropy;

        public static PascalDistribution? Fit(IEnumerable<int> samples) {
            if (!samples.Any() || samples.Any(n => n < 0)) {
                return null;
            }

            (ddouble mean, ddouble variance) = samples.Select(n => (ddouble)n).MeanVariance();
            int N = (int)Clamp(Round(mean * mean / (variance + mean)), 1, int.MaxValue);
            ddouble p = N / mean;

            try {
                return new PascalDistribution(N, p);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string Formula => "f(k; n, p) := negative_binomial(k - n; n, p)";

        public override string ToString() {
            return $"{typeof(PascalDistribution).Name}[n={N},p={P}]";
        }

    }
}
