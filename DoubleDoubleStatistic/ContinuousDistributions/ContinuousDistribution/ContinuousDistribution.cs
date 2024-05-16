using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    public abstract class ContinuousDistribution {
        public abstract ddouble PDF(ddouble x);
        public virtual ddouble CDF(ddouble x, Interval interval = Interval.Lower) => throw new NotImplementedException();
        public virtual ddouble Quantile(ddouble p, Interval interval = Interval.Lower) => throw new NotImplementedException();

        public virtual double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();
            double v = (double)Quantile(u);

            return v;
        }

        public IEnumerable<double> Sample(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return Sample(random);
                count--;
            }
        }

        public virtual bool AdditiveClosed => false;
        public virtual bool SubtractiveClosed => false;
        public virtual bool MultiplyClosed => false;
        public virtual bool Scalable => false;
        public virtual bool Shiftable => false;
        public virtual bool Symmetric => false;

        public virtual (ddouble min, ddouble max) Support =>
            (ddouble.NegativeInfinity, ddouble.PositiveInfinity);

        public virtual ddouble Mean => throw new NotImplementedException();
        public virtual ddouble Median => throw new NotImplementedException();
        public virtual ddouble Mode => throw new NotImplementedException();
        public virtual ddouble Variance => throw new NotImplementedException();
        public virtual ddouble Skewness => throw new NotImplementedException();
        public virtual ddouble Kurtosis => throw new NotImplementedException();
        public virtual ddouble Entropy => throw new NotImplementedException();
        public virtual string Formula => throw new NotFiniteNumberException();

        protected static bool InRangeUnit(ddouble v) {
            return v >= 0d && v <= 1d;
        }
    }
}