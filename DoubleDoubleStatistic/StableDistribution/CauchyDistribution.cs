using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class CauchyDistribution : StableDistribution<CauchyDistribution>,
        IAdditionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        ISubtractionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        IAdditionOperators<CauchyDistribution, ddouble, CauchyDistribution>,
        ISubtractionOperators<CauchyDistribution, ddouble, CauchyDistribution>,
        IMultiplyOperators<CauchyDistribution, ddouble, CauchyDistribution> {

        public override ddouble Mu { get; }
        public ddouble Gamma { get; }

        private readonly ddouble pdf_norm, gamma_inv, gamma_sq;

        public CauchyDistribution() : this(mu: 0, gamma: 1) { }

        public CauchyDistribution(ddouble mu, ddouble gamma) {
            ValidateLocation(mu);
            ValidateScale(gamma);

            Mu = mu;
            Gamma = gamma;

            pdf_norm = RcpPI * Gamma;
            gamma_inv = 1d / gamma;
            gamma_sq = gamma * gamma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = pdf_norm / (Square(x - Mu) + gamma_sq);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = Atan((x - Mu) * gamma_inv) * RcpPI + Point5;

                return cdf;
            }
            else {
                ddouble cdf = Atan((Mu - x) * gamma_inv) * RcpPI + Point5;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Mu + Gamma * TanPI(p - Point5);

                return x;
            }
            else {
                ddouble x = Mu - Gamma * TanPI(p - Point5);

                return x;
            }
        }

        public override bool Symmetric => true;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => Log(4 * PI * Gamma);

        public override ddouble Alpha => 1d;

        public override ddouble Beta => 0d;

        public override ddouble C => Gamma;

        public static CauchyDistribution operator +(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator -(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator +(CauchyDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Gamma);
        }

        public static CauchyDistribution operator -(CauchyDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Gamma);
        }

        public static CauchyDistribution operator *(CauchyDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Gamma * k);
        }

        public override string ToString() {
            return $"{typeof(CauchyDistribution).Name}[mu={Mu},gamma={Gamma}]";
        }

        public override string Formula => "p(x; mu, gamma) := 1 / (1 + u^2) / (pi * gamma), u = (x - mu) / gamma";
    }
}
