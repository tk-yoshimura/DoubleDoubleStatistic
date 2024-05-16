using DoubleDouble;

namespace DoubleDoubleStatistic.DirectionalDistributions {
    public abstract class DirectionalDistribution<DD, D> {
        public abstract ddouble PDF(DD v);

        public abstract D Sample(Random random);

        public IEnumerable<D> Sample(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return Sample(random);
                count--;
            }
        }

        public virtual DD Mean => throw new NotImplementedException();
        public virtual DD Mode => throw new NotImplementedException();
        public virtual DD Variance => throw new NotImplementedException();
        public virtual DD Skewness => throw new NotImplementedException();
        public virtual DD Kurtosis => throw new NotImplementedException();
        public virtual ddouble Entropy => throw new NotImplementedException();
        public virtual string Formula => throw new NotFiniteNumberException();
    }
}
