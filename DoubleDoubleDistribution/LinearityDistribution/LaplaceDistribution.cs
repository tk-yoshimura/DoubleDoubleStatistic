using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class LaplaceDistribution : LinearityDistribution<LaplaceDistribution>,
        IAdditionOperators<LaplaceDistribution, ddouble, LaplaceDistribution>,
        IMultiplyOperators<LaplaceDistribution, ddouble, LaplaceDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public LaplaceDistribution() : this(mu: 0d, sigma: 1d) { }

        public LaplaceDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;
            sigma_inv = 1d / Sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = Abs(x - Mu);

            ddouble pdf = Exp(-u * sigma_inv) * sigma_inv / 2;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x - Mu;

            if (interval == Interval.Lower) {
                ddouble cdf = IsNegative(u) ? Exp(u * sigma_inv) / 2 : 1d - Exp(-u * sigma_inv) / 2;

                return cdf;
            }
            else {
                ddouble cdf = IsPositive(u) ? Exp(-u * sigma_inv) / 2 : 1d - Exp(u * sigma_inv) / 2;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble u = 2 * p;

            if (interval == Interval.Lower) {
                ddouble x = p < 0.5d
                    ? Mu + Log(u) * Sigma
                    : Mu - Log(2d - u) * Sigma;

                return x;
            }
            else {
                ddouble x = p < 0.5d
                    ? Mu - Log(u) * Sigma
                    : Mu + Log(2d - u) * Sigma;

                return x;
            }
        }

        public override bool Symmetric => true;


        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance => 2 * Sigma * Sigma;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => 3d;

        public override ddouble Entropy => Log(2 * Sigma * E);

        public static LaplaceDistribution operator +(LaplaceDistribution dist, ddouble s) {
            return new(s + dist.Mu, dist.Sigma);
        }

        public static LaplaceDistribution operator *(LaplaceDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(LaplaceDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }
    }
}
