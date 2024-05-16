using DoubleDouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    public abstract class DiscreteDistribution {
        public abstract ddouble PMF(int k);

        public virtual int Sample(Random random) => throw new NotImplementedException();

        public IEnumerable<int> Sample(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return Sample(random);
                count--;
            }
        }

        public virtual (int min, int max) Support =>
            (0, int.MaxValue);

        public virtual ddouble Mean => throw new NotImplementedException();
        public virtual ddouble Variance => throw new NotImplementedException();
        public virtual ddouble Skewness => throw new NotImplementedException();
        public virtual ddouble Kurtosis => throw new NotImplementedException();
        public virtual ddouble Entropy => throw new NotImplementedException();
        public virtual string Formula => throw new NotFiniteNumberException();
    }
}