using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class BirnbaumSaundersDistribution : ScalableDistribution<BirnbaumSaundersDistribution>,
        IMultiplyOperators<BirnbaumSaundersDistribution, ddouble, BirnbaumSaundersDistribution> {

        public ddouble Theta { get; }
        public ddouble Alpha { get; }

        private readonly ddouble theta_inv, pdf_norm, alpha_sq;

        public BirnbaumSaundersDistribution(ddouble alpha, ddouble theta) {
            ValidateScale(theta);

            Alpha = alpha;
            Theta = theta;

            theta_inv = 1d / theta;
            alpha_sq = alpha * alpha;
            pdf_norm = 1d / (2 * Sqrt2 * Sqrt(PI / theta) * alpha);
        }

        public override ddouble PDF(ddouble x) {
            if (x <= 0d) {
                return 0d;
            }

            ddouble u = x * theta_inv, up1 = u + 1;

            ddouble pdf = pdf_norm * up1 * Exp(-up1 * up1 / (2 * alpha_sq * u));

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x)) {
                    return 1d;
                }

                ddouble cdf = Erfc((1d - u) / (Sqrt2 * Alpha * Sqrt(u))) / 2;

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble cdf = Erfc((u - 1d) / (Sqrt2 * Alpha * Sqrt(u))) / 2;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p == 0d) {
                return interval == Interval.Lower ? 0 : PositiveInfinity;
            }
            if (p == 1d) {
                return interval == Interval.Lower ? PositiveInfinity : 0;
            }

            ddouble w = Alpha * InverseErfc(p * 2), w2 = w * w;

            ddouble u;

            if (interval == Interval.Lower) {
                u = w * (w - Sqrt(w2 + 2d)) + 1d;
            }
            else {
                u = 1d / (w * (w - Sqrt(w2 + 2d)) + 1d);
            }

            ddouble x = u * Theta;

            return x;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Theta * (alpha_sq + 2d) / 2d;

        public override ddouble Median => Theta;

        public override ddouble Variance => Theta * alpha_sq * (5 * alpha_sq + 4) / 2d;

        public override ddouble Skewness => 4 * Alpha * (11 * alpha_sq + 6d) / Cube(Sqrt(5 * alpha_sq + 4d));

        public override ddouble Kurtosis => (48d + alpha_sq * (360d + alpha_sq * 633d)) / Square(5 * alpha_sq + 4d);

        public static BirnbaumSaundersDistribution operator *(BirnbaumSaundersDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Theta * k);
        }

        public override string ToString() {
            return $"{typeof(BirnbaumSaundersDistribution).Name}[alpha={Alpha},theta={Theta}]";
        }
    }
}
