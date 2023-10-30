using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class LogNormalDistribution : Distribution,
        IAdditionOperators<LogNormalDistribution, LogNormalDistribution, LogNormalDistribution>,
        ISubtractionOperators<LogNormalDistribution, LogNormalDistribution, LogNormalDistribution>,
        IMultiplyOperators<LogNormalDistribution, ddouble, LogNormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        public LogNormalDistribution() : this(mu: 0, sigma: 1) { }

        public LogNormalDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            this.Mu = mu;
            this.Sigma = sigma;

            this.sigma_sq = sigma * sigma;
            this.pdf_norm = 1d / (sigma * Sqrt(2 * PI));
            this.exp_scale = -1d / (2 * sigma_sq);
            this.erf_scale = 1d / (Sqrt2 * sigma);
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = pdf_norm * Exp(Square(x - Mu) * exp_scale);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble cdf = Ldexp(1d + Erf((x - Mu) * erf_scale), -1);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            ddouble quantile = Mu + Sigma * Sqrt2 * InverseErfc(Ldexp(p, 1));

            return quantile;
        }

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;
        public override bool Scalable => true;
        public override bool Symmetric => true;

        public override ddouble Mean => Mu;
        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;
        public override ddouble Variance => sigma_sq;
        public override ddouble Skewness => 0;
        public override ddouble Kurtosis => 0;

        public override ddouble Entropy => Log(Sigma * Sqrt(2 * PI * E));

        public static LogNormalDistribution operator +(LogNormalDistribution dist1, LogNormalDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static LogNormalDistribution operator -(LogNormalDistribution dist1, LogNormalDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static LogNormalDistribution operator *(LogNormalDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(LogNormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }
    }
}
