using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class CauchyDistribution : Distribution,
        IAdditionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        ISubtractionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution> {

        public ddouble Mu { get; }
        public ddouble Gamma { get; }

        private readonly ddouble inv_gamma, norm;

        public CauchyDistribution() : this(0, 1) { }

        public CauchyDistribution(ddouble mu, ddouble gamma) {
            ValidateLocation(mu);
            ValidateScale(gamma);

            this.Mu = mu;
            this.Gamma = gamma;

            this.inv_gamma = 1d / gamma;
            this.norm = PI * gamma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = 1d / (norm * (1d + Square((x - Mu) * inv_gamma)));

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            ddouble cdf = Atan((x - Mu) * inv_gamma) * RcpPI + Point5;

            return cdf;
        }

        public override ddouble Quantile(ddouble p) {
            ValidateProb(p);

            ddouble quantile = Mu + Gamma * TanPI(p - Point5);

            return quantile;
        }

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;

        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;

        public override ddouble Entropy => Log(4 * PI * Gamma);

        public static CauchyDistribution operator +(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator -(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public override string ToString() {
            return $"{typeof(CauchyDistribution).Name}[mu={Mu},gamma={Gamma}]";
        }
    }
}
