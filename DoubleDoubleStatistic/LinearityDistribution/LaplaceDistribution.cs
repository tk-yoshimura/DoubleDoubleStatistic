using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class LaplaceDistribution : LinearityDistribution<LaplaceDistribution>,
        IAdditionOperators<LaplaceDistribution, ddouble, LaplaceDistribution>,
        ISubtractionOperators<LaplaceDistribution, ddouble, LaplaceDistribution>,
        IMultiplyOperators<LaplaceDistribution, ddouble, LaplaceDistribution>,
        IDivisionOperators<LaplaceDistribution, ddouble, LaplaceDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public LaplaceDistribution() : this(mu: 0d, sigma: 1d) { }

        public LaplaceDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public LaplaceDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;
            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            ddouble pdf = Exp(-Abs(u)) * sigma_inv * 0.5d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x - Mu;

            if (interval == Interval.Lower) {
                ddouble cdf = IsNegative(u) ? Exp(u * sigma_inv) * 0.5d : 1d - Exp(-u * sigma_inv) * 0.5d;

                return cdf;
            }
            else {
                ddouble cdf = IsPositive(u) ? Exp(-u * sigma_inv) * 0.5d : 1d - Exp(u * sigma_inv) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble u = 2d * p;

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

        public override ddouble Variance => 2d * Sigma * Sigma;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => 3d;

        public override ddouble Entropy => Log(2d * Sigma * E);

        public static LaplaceDistribution operator +(LaplaceDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Sigma);
        }

        public static LaplaceDistribution operator -(LaplaceDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Sigma);
        }

        public static LaplaceDistribution operator *(LaplaceDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Sigma * k);
        }

        public static LaplaceDistribution operator /(LaplaceDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(LaplaceDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-|u|) / (sigma * 2), u = (x - mu) / sigma";
    }
}
