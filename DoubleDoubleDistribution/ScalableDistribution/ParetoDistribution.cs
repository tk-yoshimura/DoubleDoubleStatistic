using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ParetoDistribution : ScalableDistribution<ParetoDistribution>,
        IMultiplyOperators<ParetoDistribution, ddouble, ParetoDistribution> {

        public ddouble Xm { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm;

        public ParetoDistribution(ddouble xm, ddouble alpha) {
            Xm = xm;
            Alpha = alpha;

            pdf_norm = alpha * Pow(xm, alpha);
        }

        public override ddouble PDF(ddouble x) {
            if (x < Xm) {
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
                if (x < Xm) {
                    return 0d;
                }

                ddouble cdf = 1d - Pow(Xm / x, Alpha);

                return cdf;
            }
            else {
                if (x < Xm) {
                    return 1d;
                }

                ddouble cdf = Pow(Xm / x, Alpha);

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

            ddouble x = Xm / Pow(p, 1d / Alpha);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (Xm, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d) ? (Alpha * Xm) / (Alpha - 1d) : PositiveInfinity;

        public override ddouble Median => Xm * Pow(2d, 1d / Alpha);

        public override ddouble Mode => Xm;

        public override ddouble Variance => (Alpha > 2d)
            ? (Alpha * Xm * Xm) / (Square(Alpha - 1d) * (Alpha - 2d))
            : PositiveInfinity;

        public override ddouble Skewness => (Alpha > 3d)
            ? 2 * (Alpha + 1d) / (Alpha - 3d) * Sqrt((Alpha - 2d) / Alpha)
            : NaN;

        public override ddouble Kurtosis => (Alpha > 4d)
            ? 6d * (-2d + Alpha * (-6 + Alpha * (1d + Alpha))) / (Alpha * (Alpha - 3d) * (Alpha - 4d))
            : NaN;

        public override ddouble Entropy => 1d + Log(Xm / Alpha) + 1d / Alpha;

        public static ParetoDistribution operator *(ParetoDistribution dist, ddouble k) {
            return new(k * dist.Xm, dist.Alpha);
        }

        public override string ToString() {
            return $"{typeof(ParetoDistribution).Name}[xm={Xm},alpha={Alpha}]";
        }
    }
}
