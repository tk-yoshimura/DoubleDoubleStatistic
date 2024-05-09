using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BinaryDistribution : DiscreteDistribution {

        public BinaryDistribution() { }

        public override ddouble PMF(int k) {
            return k == 0 || k == 1 ? 0.5d : 0d;
        }

        public override int Sample(Random random) {
            return (int)unchecked(random.NextInt64() & 1L);
        }

        public override (int min, int max) Support => (0, 1);

        public override ddouble Mean => 0.5d;

        public override ddouble Variance => 0.25d;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => -2d;

        public override ddouble Entropy => Ln2;

        public override string Formula => "f(k) := 1 / 2";

        public override string ToString() {
            return $"{typeof(BinaryDistribution).Name}[]";
        }

    }
}
