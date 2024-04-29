using DoubleDouble;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;
using Complex = DoubleDoubleComplex.Complex;

namespace DoubleDoubleStatistic {
    public class VoigtDistribution : ScalableDistribution<VoigtDistribution>,
        IAdditionOperators<VoigtDistribution, VoigtDistribution, VoigtDistribution>,
        IMultiplyOperators<VoigtDistribution, ddouble, VoigtDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, z_scale, zr, cdf_limit;

        public bool EnableCDFErrorException { get; set; } = false;

        public VoigtDistribution() : this(gamma: 1, sigma: 1) { }

        private const int cache_samples = 512;
        private CDFSegmentCache cdf_cache;
        private QuantileBuilder quantile_lower_builder = null, quantile_upper_builder = null;

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

            if (IsInfinity(x)) {
                if (x > 0d) {
                    return interval == Interval.Lower ? 1d : 0d;
                }
                else { 
                    return interval == Interval.Lower ? 0d : 1d;
                }
            }

            ddouble p = PDF(x);
            if (IsZero(p)) {
                return x < 0d ? 0d : 1d;
            }

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

            this.cdf_cache ??= new CDFSegmentCache(
                0d, 1d,
                f,
                samples: cache_samples
            );

            ddouble t = 1d / (Abs(x) + 1d);

            ddouble cdf = (interval == Interval.Lower) ^ (x < 0d)
                ? Min(1d, 0.5d + cdf_cache.Upper(t) * pdf_norm)
                : Min(0.5d, cdf_cache.Lower(t) * pdf_norm);

            return cdf;
        }

        public override bool AdditiveClosed => true;

        public override bool Symmetric => true;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static VoigtDistribution operator +(VoigtDistribution dist1, VoigtDistribution dist2) {
            return new(dist1.Gamma + dist2.Gamma, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static VoigtDistribution operator *(VoigtDistribution dist, ddouble k) {
            return new(dist.Gamma * k, dist.Sigma * k);
        }

        public override string Formula => "p(x; gamma, sigma) := Re[erfcx((x + i * gamma) / (sigma * sqrt(2)))] / (sigma * sqrt(2 * pi))";

        public override string ToString() {
            return $"{typeof(VoigtDistribution).Name}[gamma={Gamma},sigma={Sigma}]";
        }
    }
}
