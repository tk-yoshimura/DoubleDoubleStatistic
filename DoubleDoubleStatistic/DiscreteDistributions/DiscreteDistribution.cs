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

        protected static void ValidateScale(ddouble scale) {
            if (!(scale > 0d && ddouble.IsFinite(scale))) {
                throw new ArgumentOutOfRangeException(nameof(scale), "Invalid scale parameter.");
            }
        }

        protected static void ValidateLocation(ddouble location) {
            if (!ddouble.IsFinite(location)) {
                throw new ArgumentOutOfRangeException(nameof(location), "Invalid location parameter.");
            }
        }

        protected static void ValidateLocation(ddouble location, Func<ddouble, bool> condition) {
            if (!ddouble.IsFinite(location) || !condition(location)) {
                throw new ArgumentOutOfRangeException(nameof(location), "Invalid location parameter.");
            }
        }

        protected static void ValidateShape(ddouble shape, Func<ddouble, bool> condition) {
            if (!condition(shape)) {
                throw new ArgumentOutOfRangeException(nameof(shape), "Invalid shape parameter.");
            }
        }

        protected static void ValidateShape(int shape, Func<int, bool> condition) {
            if (!condition(shape)) {
                throw new ArgumentOutOfRangeException(nameof(shape), "Invalid shape parameter.");
            }
        }
    }
}