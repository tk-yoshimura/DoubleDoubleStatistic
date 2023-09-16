using DoubleDouble;

namespace DoubleDoubleDistribution {
    public abstract class Distribution {
        public abstract ddouble PDF(ddouble x);
        public virtual ddouble CDF(ddouble x) => throw new NotImplementedException();
        public virtual ddouble Quantile(ddouble p) => throw new NotImplementedException();

        public virtual bool AdditiveClosed => false;

        public virtual (ddouble min, ddouble max) Domain => 
            (ddouble.NegativeInfinity, ddouble.PositiveInfinity);

        public virtual ddouble Mean => ddouble.NaN;
        public virtual ddouble Median => ddouble.NaN;
        public virtual ddouble Mode => ddouble.NaN;
        public virtual ddouble Variance => ddouble.NaN;
        public virtual ddouble Skewness => ddouble.NaN;
        public virtual ddouble Kurtosis => ddouble.NaN;

        public virtual ddouble Entropy => ddouble.NaN;

        public static void ValidateProb(ddouble p) {
            if (!(p >= 0d && p <= 1d)) {
                throw new ArgumentOutOfRangeException(nameof(p), "Invalid probably.");
            }
        }

        public static void ValidateScale(ddouble scale) {
            if (!(scale > 0d && ddouble.IsFinite(scale))) {
                throw new ArgumentOutOfRangeException(nameof(scale), "Invalid scale parameter.");
            }
        }

        public static void ValidateLocation(ddouble location) {
            if (!ddouble.IsFinite(location)) {
                throw new ArgumentOutOfRangeException(nameof(location), "Invalid location parameter.");
            }
        }
    }
}