using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class SkewNormalDistribution : LinearityDistribution<SkewNormalDistribution>,
        IAdditionOperators<SkewNormalDistribution, ddouble, SkewNormalDistribution>,
        ISubtractionOperators<SkewNormalDistribution, ddouble, SkewNormalDistribution>,
        IMultiplyOperators<SkewNormalDistribution, ddouble, SkewNormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm, sigma_inv, s;

        public SkewNormalDistribution() : this(alpha: 0, mu: 0, sigma: 1) { }

        public SkewNormalDistribution(ddouble alpha, ddouble mu, ddouble sigma) {
            ValidateShape(alpha, IsFinite);
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;
            Alpha = alpha;

            pdf_norm = 2d / (Square(sigma) * Sqrt(2d * PI));
            sigma_inv = 1d / sigma;

            s = Alpha / Hypot(1, alpha);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;

            ddouble pdf = pdf_norm * Exp(-u * u * 0.5d) * Erfc(-Alpha * u / Sqrt2);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;

            if (interval == Interval.Lower) {
                ddouble cdf = Erfc(-Alpha * u / Sqrt2) / 2 - 2 * OwenT(u, Alpha);
                cdf = IsFinite(cdf) ? Clamp(cdf, 0d, 1d) : (x < Mu) ? 0d : 1d;

                return cdf;
            }
            else {
                ddouble cdf = Erfc(Alpha * u / Sqrt2) / 2 + 2 * OwenT(u, Alpha);
                cdf = IsFinite(cdf) ? Clamp(cdf, 0d, 1d) : (x < Mu) ? 1d : 0d;

                return cdf;
            }
        }

        public override ddouble Mean =>
            Mu + Sigma * s * Sqrt(2d / PI);

        public override ddouble Variance =>
            Sigma * Sigma * (1d - 2d * s * s / PI);

        public override ddouble Skewness =>
            (4d - PI) / 2 * Cube(s * Sqrt(2 / PI)) / Cube(Sqrt(1 - 2 * s * s / PI));

        public override ddouble Kurtosis =>
            2 * (PI - 3d) * Square(Square(s * Sqrt(2 / PI))) / Square(1 - 2 * s * s / PI);

        public static SkewNormalDistribution operator +(SkewNormalDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Sigma);
        }

        public static SkewNormalDistribution operator -(SkewNormalDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Sigma);
        }

        public static SkewNormalDistribution operator *(SkewNormalDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Sigma * k);
        }

        public override string ToString() {
            return $"{typeof(SkewNormalDistribution).Name}[alpha={Alpha},mu={Mu},sigma={Sigma}]";
        }
    }
}
