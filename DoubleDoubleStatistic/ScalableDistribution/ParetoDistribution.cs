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

            pdf_norm = alpha / k;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = K / x;

            if (IsNegative(x) || u > 1d) {
                return 0d;
            }
            if (IsNaN(u)) {
                return NaN;
            }

            ddouble pdf = pdf_norm * Pow(u, Alpha + 1d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = K / x;

            if (interval == Interval.Lower) {
                if (IsNegative(x) || u > 1d) {
                    return 0d;
                }

                ddouble cdf = 1d - Pow(u, Alpha);

                return cdf;
            }
            else {
                if (IsNegative(x) || u > 1d) {
                    return 1d;
                }

                ddouble cdf = Pow(u, Alpha);

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

        public override string Formula => "p(x; k, alpha) := u^(- alpha - 1) * (alpha / k), u = x / K";
    }
}
