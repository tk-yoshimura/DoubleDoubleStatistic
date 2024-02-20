using DoubleDouble;
using DoubleDoubleIntegrate;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;
using Complex = DoubleDoubleComplex.Complex;

namespace DoubleDoubleDistribution {
    public class VoigtDistribution : ScalableDistribution<VoigtDistribution>,
        IAdditionOperators<VoigtDistribution, VoigtDistribution, VoigtDistribution>,
        IMultiplyOperators<VoigtDistribution, ddouble, VoigtDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, z_scale, zr, cdf_limit;

        public bool EnableCDFErrorException { get; set; } = false;

        public VoigtDistribution() : this(gamma: 1, sigma: 1) { }

        public VoigtDistribution(ddouble gamma, ddouble sigma) {
            ValidateScale(gamma);
            ValidateScale(sigma);

            Gamma = gamma;
            Sigma = sigma;

            pdf_norm = 1d / (sigma * Sqrt(2 * PI));
            z_scale = -1d / (Sqrt2 * sigma);
            zr = -gamma * z_scale;
            cdf_limit = gamma * sigma * Sqrt(2 * RcpPI);
        }

        public override ddouble PDF(ddouble x) {
            if (IsInfinity(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            Complex z = (zr, x * z_scale);

            ddouble pdf = Complex.Erfcx(z).R * pdf_norm;
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return CDF(-x, Interval.Lower);
            }

            ddouble p = PDF(x);
            if (IsZero(p)) {
                return x < 0d ? 0d : 1d;
            }

            ddouble eps = 1e-27 * pdf_norm;

            ddouble u = 1d / (Abs(x) + 1d);

            ddouble error, cdf;

            if (u < 0.5d) {
                ddouble f(ddouble t) {
                    if (IsZero(t)) {
                        return cdf_limit;
                    }

                    ddouble t_inv = 1d / t;
                    ddouble x = (1d - t) * t_inv;

                    Complex z = (zr, x * z_scale);
                    ddouble pdf = Complex.Erfcx(z).R;

                    ddouble y = pdf * t_inv * t_inv;

                    return y;
                };

                (ddouble value, error, long eval_points) = 
                    GaussKronrodIntegral.AdaptiveIntegrate(f, 0d, u, eps, order: GaussKronrodOrder.G31K63, maxdepth: 6);
                value = Max(0d, value) * pdf_norm;

                cdf = x < 0d ? value : 1d - value;

                Debug.WriteLine($"evals: {eval_points}");
            }
            else {
                ddouble f(ddouble t) {
                    ddouble t_inv = 1d / t;
                    ddouble x = (1d - t) * t_inv;

                    Complex z = (zr, x * z_scale);
                    ddouble pdf = Complex.Erfcx(z).R;

                    ddouble y = pdf * t_inv * t_inv;

                    return y;
                };

                (ddouble value, error, long eval_points) = 
                    GaussKronrodIntegral.AdaptiveIntegrate(f, u, 1d, eps, order: GaussKronrodOrder.G31K63, maxdepth: 6);
                value = Max(0d, value) * pdf_norm;

                cdf = x < 0d ? 0.5d - value : 0.5d + value;

                Debug.WriteLine($"evals: {eval_points}");
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

            Complex z = (zr, x * z_scale);
            ddouble pdf = Complex.Erfcx(z).R;

            ddouble y = pdf * t_inv * t_inv;

            return y;
        }

        public override bool AdditiveClosed => true;

        public override bool Symmetric => true;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public static VoigtDistribution operator +(VoigtDistribution dist1, VoigtDistribution dist2) {
            return new(dist1.Gamma + dist2.Gamma, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static VoigtDistribution operator *(VoigtDistribution dist, ddouble k) {
            return new(dist.Gamma * k, dist.Sigma * k);
        }

        public override string ToString() {
            return $"{typeof(VoigtDistribution).Name}[gamma={Gamma},sigma={Sigma}]";
        }
    }
}
