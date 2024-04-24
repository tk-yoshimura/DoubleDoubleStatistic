using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class ParetoDistribution : ScalableDistribution<ParetoDistribution>,
        IMultiplyOperators<ParetoDistribution, ddouble, ParetoDistribution> {

        public ddouble K { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm;

        public ParetoDistribution(ddouble k, ddouble alpha) {
            K = k;
            Alpha = alpha;

            pdf_norm = alpha * Pow(k, alpha);
        }

        public override ddouble PDF(ddouble x) {
            if (x < K) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble pdf = pdf_norm / Pow(x, Alpha + 1d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < K) {
                    return 0d;
                }

                ddouble cdf = 1d - Pow(K / x, Alpha);

                return cdf;
            }
            else {
                if (x < K) {
                    return 1d;
                }

                ddouble cdf = Pow(K / x, Alpha);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return Quantile(1d - p, Interval.Upper);
            }

            ddouble x = K / Pow(p, 1d / Alpha);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (K, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d) ? (Alpha * K) / (Alpha - 1d) : PositiveInfinity;

        public override ddouble Median => K * Pow(2d, 1d / Alpha);

        public override ddouble Mode => K;

        public override ddouble Variance => (Alpha > 2d)
            ? (Alpha * K * K) / (Square(Alpha - 1d) * (Alpha - 2d))
            : PositiveInfinity;

        public override ddouble Skewness => (Alpha > 3d)
            ? 2 * (Alpha + 1d) / (Alpha - 3d) * Sqrt((Alpha - 2d) / Alpha)
            : NaN;

        public override ddouble Kurtosis => (Alpha > 4d)
            ? 6d * (-2d + Alpha * (-6 + Alpha * (1d + Alpha))) / (Alpha * (Alpha - 3d) * (Alpha - 4d))
            : NaN;

        public override ddouble Entropy => 1d + Log(K / Alpha) + 1d / Alpha;

        public static ParetoDistribution operator *(ParetoDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.K * k);
        }

        public override string ToString() {
            return $"{typeof(ParetoDistribution).Name}[k={K},alpha={Alpha}]";
        }
    }
}
