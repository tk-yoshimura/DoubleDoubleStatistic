using DoubleDouble;

namespace DoubleDoubleDistribution {
    public abstract class Distribution {
        public abstract ddouble PDF(ddouble x);
        public virtual ddouble CDF(ddouble x) => throw new NotImplementedException();
        public virtual ddouble Quantile(ddouble p) => throw new NotImplementedException();

        public virtual bool AdditiveClosed => false;
        public virtual bool SubtractiveClosed => false;
        public virtual bool Scalable => false;
        public virtual bool Symmetric => false;

        public virtual (ddouble min, ddouble max) Support =>
            (ddouble.NegativeInfinity, ddouble.PositiveInfinity);

        public virtual ddouble Mean => ddouble.NaN;
        public virtual ddouble Median => ddouble.NaN;
        public virtual ddouble Mode => ddouble.NaN;
        public virtual ddouble Variance => ddouble.NaN;
        public virtual ddouble Skewness => ddouble.NaN;
        public virtual ddouble Kurtosis => ddouble.NaN;

        public virtual ddouble Entropy => ddouble.NaN;

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