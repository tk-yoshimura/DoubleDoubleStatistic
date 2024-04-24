using DoubleDouble;

namespace DoubleDoubleDistribution {
    public abstract class ContinuousDistribution {
        public abstract ddouble PDF(ddouble x);
        public virtual ddouble CDF(ddouble x, Interval interval = Interval.Lower) => throw new NotImplementedException();
        public virtual ddouble Quantile(ddouble p, Interval interval = Interval.Lower) => throw new NotImplementedException();

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

        protected static bool InRangeUnit(ddouble v) {
            return v >= 0d && v <= 1d;
        }

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