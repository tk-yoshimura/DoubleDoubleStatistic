using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class GaussKuzminDistribution : DiscreteDistribution {
        const int random_gen_max_index = 65536;

        private Roulette roulette = null;

        public GaussKuzminDistribution() {
        }

        public override ddouble PMF(int k) {
            return k >= 1 ? -Log2(1d - Square((ddouble)1d / (k + 1d))) : 0d;
        }

        public override int Sample(Random random) {
            if (roulette is null) {
                SetupRoulette();
            }

            return roulette.NextIndex(random);
        }

        public void SetupRoulette() {
            List<double> probs = [];
            ddouble ps = 0d;
            for (int k = 0; k <= random_gen_max_index && (double)ps < 1d; k++) {
                ddouble p = PMF(k);
                ps += p;

                probs.Add((double)p);
            }

            roulette = new(probs);
        }

        public override (int min, int max) Support => (1, int.MaxValue);

        public override ddouble Mean => PositiveInfinity;

        public override ddouble Variance => PositiveInfinity;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => "2.379246769061239570532";

        public override string Formula => "f(k) := -log2(1 - 1 / (k + 1)^2)";

        public override string ToString() {
            return $"{typeof(GaussKuzminDistribution).Name}[]";
        }

    }
}
