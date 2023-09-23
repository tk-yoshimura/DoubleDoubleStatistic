using DoubleDouble;
using DoubleDoubleComplex;
using DoubleDoubleIntegrate;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class VoigtDistribution : Distribution,
        System.Numerics.IAdditionOperators<VoigtDistribution, VoigtDistribution, VoigtDistribution>,
        System.Numerics.IMultiplyOperators<VoigtDistribution, ddouble, VoigtDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private readonly ddouble norm, inv_scale, i, cdf_limit;

        public bool EnableCDFErrorException { get; set; } = false;

        public VoigtDistribution() : this(gamma: 1, sigma: 1) { }

        public VoigtDistribution(ddouble gamma, ddouble sigma) {
            ValidateScale(gamma);
            ValidateScale(sigma);

            this.Gamma = gamma;
            this.Sigma = sigma;

            this.norm = 1d / (sigma * Sqrt(2 * PI));
            this.inv_scale = -1d / (Sqrt2 * sigma);
            this.i = -gamma * inv_scale;
            this.cdf_limit = gamma * sigma * Sqrt(2 * RcpPI);
        }

        public override ddouble PDF(ddouble x) {
            Complex z = (i, x * inv_scale);

            ddouble pdf = Complex.Erfcx(z).R * norm;
            pdf = IsNaN(pdf) ? 0d : pdf;

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            ddouble p = PDF(x);
            if (IsZero(p)) {
                return x < 0d ? 0d : 1d;
            }

            ddouble eps = 1e-27 * norm;

            ddouble u = 1d / (Abs(x) + 1d);

            ddouble error, cdf;

            if (u < 0.5d) {
                ddouble f(ddouble t) {
                    if (IsZero(t)) {
                        return cdf_limit;
                    }

                    ddouble t_inv = 1d / t;
                    ddouble x = (1d - t) * t_inv;

                    Complex z = (i, x * inv_scale);
                    ddouble pdf = Complex.Erfcx(z).R;

                    ddouble y = pdf * t_inv * t_inv;

                    return y;
                };

                (ddouble value, error) = GaussKronrodIntegral.AdaptiveIntegrate(f, Zero, u, eps, depth: 12);
                value = Max(0d, value) * norm;

                cdf = x < 0d ? value : 1d - value;
            }
            else {
                ddouble f(ddouble t) {
                    ddouble t_inv = 1d / t;
                    ddouble x = (1d - t) * t_inv;

                    Complex z = (i, x * inv_scale);
                    ddouble pdf = Complex.Erfcx(z).R;

                    ddouble y = pdf * t_inv * t_inv;

                    return y;
                };

                (ddouble value, error) = GaussKronrodIntegral.AdaptiveIntegrate(f, u, One, eps, depth: 12);
                value = Max(0d, value) * norm;

                cdf = x < 0d ? 0.5d - value : 0.5d + value;
            }

            if (EnableCDFErrorException && !(error < eps)) {
                throw new ArithmeticException("CDF integrate not convergence.");
            }

            return cdf;
        }

        internal ddouble Integrand(ddouble t) {
            if (IsZero(t)) {
                return cdf_limit;
            }

            ddouble t_inv = 1d / t;
            ddouble x = (1d - t) * t_inv;

            Complex z = (i, x * inv_scale);
            ddouble pdf = Complex.Erfcx(z).R;

            ddouble y = pdf * t_inv * t_inv;

            return y;
        }

        public override bool AdditiveClosed => true;
        public override bool Scalable => true;
        public override bool Symmetric => true;

        public override ddouble Median => 0;
        public override ddouble Mode => 0;

        public static VoigtDistribution operator +(VoigtDistribution dist1, VoigtDistribution dist2) {
            return new(dist1.Gamma + dist2.Gamma, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static VoigtDistribution operator *(VoigtDistribution dist, ddouble k) {
            return new(k * dist.Gamma, k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(VoigtDistribution).Name}[gamma={Gamma},sigma={Sigma}]";
        }
    }
}
