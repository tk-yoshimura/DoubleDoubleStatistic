using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BernoulliDistribution : DiscreteDistribution {

        public ddouble P { get; }

        private readonly ddouble p, q;

        public BernoulliDistribution(ddouble p) {
            ValidateShape(p, p => p >= 0d && p <= 1d);

            P = p;
            this.p = p;
            q = 1d - p;
        }

        public override ddouble PMF(int k) {
            return k == 0 ? q : k == 1 ? p : 0d;
        }

        public override int Sample(Random random) {
            return random.NextUniformOpenInterval1() < p ? 1 : 0;
        }

        public override (int min, int max) Support => (0, 1);

        public override ddouble Mean => p;

        public override ddouble Variance => p * q;

        public override ddouble Skewness => (q - p) / Sqrt(p * q);

        public override ddouble Kurtosis => 1d / (p * q) - 6d;

        public override ddouble Entropy => -(q * Log(q) + p * Log(p));

        public override string Formula => "f(k; p) := if k = 0 then 1 - p else p";

        public override string ToString() {
            return $"{typeof(BernoulliDistribution).Name}[c={P}]";
        }

    }
}
