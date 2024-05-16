using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class CategoricalDistribution : DiscreteDistribution {

        public ReadOnlyCollection<ddouble> Probs { get; }

        public int N => Probs.Count;

        private readonly Roulette roulette;

        public CategoricalDistribution(params ddouble[] probs) : this(new ReadOnlyCollection<ddouble>(probs)) { }

        public CategoricalDistribution(IEnumerable<ddouble> probs) : this(new ReadOnlyCollection<ddouble>(probs.ToArray())) { }

        public CategoricalDistribution(ReadOnlyCollection<ddouble> probs) {
            ddouble prob_sum = probs.Sum();

            ParamAssert.ValidateShape(nameof(probs), probs.All(p => p >= 0d && IsFinite(p)) && IsFinite(prob_sum));

            Probs = new ReadOnlyCollection<ddouble>(probs.Select(p => p / prob_sum).ToArray());

            roulette = new(probs.Select(p => (double)p));
        }

        public override ddouble PMF(int k) {
            return k >= 0 && k < N ? Probs[k] : 0d;
        }

        public override int Sample(Random random) {
            return roulette.NextIndex(random);
        }

        public override (int min, int max) Support => (0, N - 1);

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            Probs.Select((p, i) => p * i).Sum();

        private ddouble? variance = null;
        public override ddouble Variance {
            get {
                if (variance is not null) {
                    return variance.Value;
                }

                ddouble mean = Mean;

                return variance ??= Probs.Select((p, i) => p * Square(i - mean)).Sum();
            }
        }

        private ddouble? skewness = null;
        public override ddouble Skewness {
            get {
                if (skewness is not null) {
                    return skewness.Value;
                }

                ddouble mean = Mean, variance = Variance;

                return skewness ??= Probs.Select((p, i) => p * Cube(i - mean)).Sum() / ExMath.Pow3d2(variance);
            }
        }

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis {
            get {
                if (kurtosis is not null) {
                    return kurtosis.Value;
                }

                ddouble mean = Mean, variance = Variance;

                return kurtosis ??= Probs.Select((p, i) => p * Square(Square(i - mean))).Sum() / Square(variance) - 3d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, N - 1);

        public override string Formula => "f(k, P) := p_i";

        public override string ToString() {
            return $"{typeof(CategoricalDistribution).Name}[P]";
        }

    }
}
