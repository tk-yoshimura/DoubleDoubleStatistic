using DoubleDouble;
using DoubleDoubleComplex;
using DoubleDoubleIntegrate;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class VoigtDistribution : Distribution,
        System.Numerics.IAdditionOperators<VoigtDistribution, VoigtDistribution, VoigtDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private readonly ddouble norm, inv_scale, i, peak;
        private readonly Distribution? limit_distribution = null;

        public VoigtDistribution() : this(1, 1) { }

        public VoigtDistribution(ddouble gamma, ddouble sigma) {
            ValidateScale(gamma);
            ValidateScale(sigma);

            this.Gamma = gamma;
            this.Sigma = sigma;

            this.norm = 1d / (sigma * Sqrt(2 * PI));
            this.inv_scale = -1d / (Sqrt2 * sigma);
            this.i = -gamma * inv_scale;

            if (Gamma == 0d) {
                limit_distribution = new NormalDistribution(mu: 0, sigma: sigma);
            }
            else if (Sigma == 0d) { 
                limit_distribution = new CauchyDistribution(mu: 0, gamma: gamma);
            }

            this.peak = PDF(0);
        }

        public override ddouble PDF(ddouble x) {
            if (limit_distribution is not null) {
                return limit_distribution.PDF(x);
            }

            Complex z = (i, x * inv_scale);

            ddouble pdf = Complex.Erfcx(z).R * norm;

            return IsNaN(pdf) ? 0d : pdf;
        }

        public override ddouble CDF(ddouble x) {
            if (limit_distribution is not null) {
                return limit_distribution.CDF(x);
            }

            ddouble p = PDF(x);
            if (IsZero(p)) {
                return x < 0d ? 0d : 1d;
            }

            const double eps = 1e-27;

            if ((double)p > (double)peak * 0.125) {
                if (x < 0d) {
                    (ddouble value, ddouble error) = GaussKronrodIntegral.AdaptiveIntegrate(PDF, x, 0d, eps, depth: 10);
                    ddouble cdf = 0.5d - value;

                    return cdf;
                }
                else {
                    (ddouble value, ddouble error) = GaussKronrodIntegral.AdaptiveIntegrate(PDF, 0d, x, eps, depth: 10);
                    ddouble cdf = 0.5d + value;

                    return cdf;
                }
            }
            else { 
                if (x < 0d) {
                    (ddouble value, ddouble error) = GaussKronrodIntegral.AdaptiveIntegrate(PDF, NegativeInfinity, x, eps, depth: 10);
                    ddouble cdf = value;

                    return cdf;
                }
                else {
                    (ddouble value, ddouble error) = GaussKronrodIntegral.AdaptiveIntegrate(PDF, x, PositiveInfinity, eps, depth: 10);
                    ddouble cdf = 1d - value;

                    return cdf;
                }
            }
        }

        public override bool AdditiveClosed => true;

        public override ddouble Median => 0;
        public override ddouble Mode => 0;

        public static VoigtDistribution operator +(VoigtDistribution dist1, VoigtDistribution dist2) {
            return new(dist1.Gamma + dist2.Gamma, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public override string ToString() {
            return $"{typeof(VoigtDistribution).Name}[gamma={Gamma},sigma={Sigma}]";
        }
    }
}
