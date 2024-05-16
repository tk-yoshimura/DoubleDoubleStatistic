using Algebra;
using DoubleDouble;

namespace DoubleDoubleStatistic.MultiVariateDistributions {
    public abstract class MultiVariateDistribution {
        public abstract ddouble PDF(Vector v);

        public abstract double[] Sample(Random random);

        public IEnumerable<double[]> Sample(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return Sample(random);
                count--;
            }
        }

        public virtual Vector Mean => throw new NotImplementedException();
        public virtual Vector Mode => throw new NotImplementedException();
        public virtual Matrix Covariance => throw new NotImplementedException();
        public virtual ddouble Entropy => throw new NotImplementedException();
        public virtual string Formula => throw new NotFiniteNumberException();
    }
}
