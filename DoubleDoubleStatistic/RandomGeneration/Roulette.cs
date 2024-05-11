using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Collections.ObjectModel;

namespace DoubleDoubleStatistic.RandomGeneration {
    public class Roulette {
        private readonly double prob_sum;
        private readonly ReadOnlyCollection<double> table;

        public Roulette(IEnumerable<double> probs) : this(new ReadOnlyCollection<double>(probs.ToArray())) { }

        public Roulette(ReadOnlyCollection<double> probs) {
            if (probs.Count < 1 || !probs.All(p => p >= 0d && double.IsFinite(p)) || !double.IsFinite(probs.Sum())) {
                throw new ArgumentException($"Invalid array: '{nameof(probs)}'", nameof(probs));
            }

            List<double> table = [0];

            ddouble s = 0d;
            foreach (double prob in probs) {
                s += prob;
                table.Add((double)s);
            }

            prob_sum = (double)s;
            this.table = new(table);
        }

        public int NextIndex(Random random) {
            double p = random.NextUniformOpenInterval1() * prob_sum;

            int index = Indexer.BisectionSearch(p, table);

            return index;
        }

        public IEnumerable<int> NextIndex(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return NextIndex(random);
                count--;
            }
        }
    }
}
