using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class CosineDistribution : LinearityDistribution<CosineDistribution>,
        IAdditionOperators<CosineDistribution, ddouble, CosineDistribution>,
        IMultiplyOperators<CosineDistribution, ddouble, CosineDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public CosineDistribution() : this(mu: 0, sigma: 1) { }

        public CosineDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;

            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            if (Abs(u) >= 1d) {
                return 0d;
            }

            ddouble pdf = (1 + CosPI(u)) / 2;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (interval == Interval.Lower) {
                if (u <= -1d) {
                    return 0d;
                }
                if (u >= 1d) {
                    return 1d;
                }

                ddouble cdf = (1 + u + SinPI(u) / PI) / 2;

                return cdf;
            }
            else {
                if (u <= -1d) {
                    return 1d;
                }
                if (u >= 1d) {
                    return 0d;
                }

                ddouble cdf = (1 - u - SinPI(u) / PI) / 2;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p < 1e-28) {
                return (interval == Interval.Lower) ? Mu - Sigma : Mu + Sigma;
            }
            if ((1d - p) < 1e-28) {
                return (interval == Interval.Lower) ? Mu + Sigma : Mu - Sigma;
            }
            if (Abs(p - 0.5d) < 1e-28) {
                return Mu;
            }

            ddouble u = KeplerE(PI * p * 2, 1d) * RcpPI - 1d;

            ddouble x = (interval == Interval.Lower) ? (Mu + Sigma * u) : (Mu - Sigma * u);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (Mu - Sigma, Mu + Sigma);

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance =>
            Square(Sigma) * (1d - 6d / (PI * PI)) / 3d;
        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis =>
            6d * (90d - Square(Square(PI))) / (5d * Square(PI * PI - 6d));

        public override ddouble Entropy => -0.5d + 2 * Ln2 + Log(Sigma);

        public static CosineDistribution operator +(CosineDistribution dist, ddouble s) {
            return new(s + dist.Mu, dist.Sigma);
        }

        public static CosineDistribution operator *(CosineDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(CosineDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }
    }
}
