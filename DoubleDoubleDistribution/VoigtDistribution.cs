using DoubleDouble;
using DoubleDoubleComplex;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class VoigtDistribution : Distribution, 
        System.Numerics.IAdditionOperators<VoigtDistribution, VoigtDistribution, VoigtDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private readonly ddouble norm, inv_scale, i;

        public VoigtDistribution() : this(1, 1) { }

        public VoigtDistribution(ddouble gamma, ddouble sigma) {
            ValidateScale(gamma);
            ValidateScale(sigma);

            this.Gamma = gamma;
            this.Sigma = sigma;

            this.norm = 1d / (sigma * Sqrt(2 * PI));
            this.inv_scale = -1d / (Sqrt2 * sigma);
            this.i = -gamma * inv_scale;
        }

        public override ddouble PDF(ddouble x) {
            Complex z = (i, x * inv_scale);

            ddouble pdf = Complex.Erfcx(z).R * norm;

            return pdf;
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
