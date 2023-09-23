using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class CauchyDistribution : Distribution,
        IAdditionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        ISubtractionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        IMultiplyOperators<CauchyDistribution, ddouble, CauchyDistribution> {

        public ddouble Mu { get; }
        public ddouble Gamma { get; }

        private readonly ddouble gamma_inv, gamma_sq;

        public CauchyDistribution() : this(mu: 0, gamma: 1) { }

        public CauchyDistribution(ddouble mu, ddouble gamma) {
            ValidateLocation(mu);
            ValidateScale(gamma);

            this.Mu = mu;
            this.Gamma = gamma;

            this.gamma_inv = 1d / gamma;
            this.gamma_sq = gamma * gamma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = RcpPI * Gamma / (Square(x - Mu) + gamma_sq);

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            ddouble cdf = Atan((x - Mu) * gamma_inv) * RcpPI + Point5;

            return cdf;
        }

        public override ddouble Quantile(ddouble p) {
            ValidateProb(p);

            ddouble quantile = Mu + Gamma * TanPI(p - Point5);

            return quantile;
        }

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;
        public override bool Scalable => true;

        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;

        public override ddouble Entropy => Log(4 * PI * Gamma);

        public static CauchyDistribution operator +(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator -(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator *(CauchyDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.Gamma);
        }

        public override string ToString() {
            return $"{typeof(CauchyDistribution).Name}[mu={Mu},gamma={Gamma}]";
        }
    }
}
