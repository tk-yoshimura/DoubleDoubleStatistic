using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class NormalDistribution : Distribution,
        IAdditionOperators<NormalDistribution, NormalDistribution, NormalDistribution>,
        ISubtractionOperators<NormalDistribution, NormalDistribution, NormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_sq, norm, exp_scale, erf_scale;

        public NormalDistribution() : this(0, 1) { }

        public NormalDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            this.Mu = mu;
            this.Sigma = sigma;

            this.sigma_sq = sigma * sigma;
            this.norm = 1d / (sigma * Sqrt(2 * PI));
            this.exp_scale = -1d / (2 * sigma_sq);
            this.erf_scale = 1d / (Sqrt2 * sigma);
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = norm * Exp(Square(x - Mu) * exp_scale);

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            ddouble cdf = Ldexp(1d + Erf((x - Mu) * erf_scale), -1);

            return cdf;
        }

        public override ddouble Quantile(ddouble p) {
            ValidateProb(p);

            ddouble quantile = Mu + Sigma * Sqrt2 * InverseErfc(Ldexp(p, 1));

            return quantile;
        }

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;

        public override ddouble Mean => Mu;
        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;
        public override ddouble Variance => sigma_sq;
        public override ddouble Skewness => 0;
        public override ddouble Kurtosis => 0;

        public override ddouble Entropy => Log(Sigma * Sqrt(2 * PI * E));

        public static NormalDistribution operator +(NormalDistribution dist1, NormalDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static NormalDistribution operator -(NormalDistribution dist1, NormalDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public override string ToString() {
            return $"{typeof(NormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }
    }
}
