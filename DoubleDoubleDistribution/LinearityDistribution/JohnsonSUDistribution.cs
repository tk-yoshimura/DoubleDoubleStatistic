using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class JohnsonSUDistribution : LinearityDistribution<JohnsonSUDistribution>,
        IMultiplyOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution>,
        IAdditionOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution>,
        ISubtractionOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution> {

        public ddouble Gamma { get; }
        public ddouble Delta { get; }
        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_inv;

        public JohnsonSUDistribution(ddouble gamma, ddouble delta, ddouble mu, ddouble sigma) {
            ValidateShape(gamma, IsFinite);
            ValidateScale(delta);
            ValidateLocation(mu);
            ValidateScale(sigma);

            Gamma = gamma;
            Delta = delta;
            Mu = mu;
            Sigma = sigma;

            pdf_norm = Delta / (Sigma * Sqrt(2 * PI));
            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;

            ddouble pdf = pdf_norm / Hypot(1, u) * Exp(-0.5d * Square(Gamma + Delta * Arsinh(u)));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;
            ddouble v = Gamma + Delta * Arsinh(u);

            if (interval == Interval.Lower) {
                ddouble cdf = Erfc(-v / Sqrt2) / 2;

                return cdf;
            }
            else {
                ddouble cdf = Erfc(v / Sqrt2) / 2;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble v = interval == Interval.Lower
                ? -InverseErfc(p * 2) * Sqrt2
                : InverseErfc(p * 2) * Sqrt2;

            ddouble u = Sinh((v - Gamma) / Delta);
            ddouble x = Mu + Sigma * u;

            return x;
        }

        public override ddouble Mean =>
            Mu - Sigma * Exp(1d / (2 * Square(Delta))) * Sinh(Gamma / Delta);

        public override ddouble Median =>
            Mu + Sigma * Sinh(-Gamma / Delta);

        public override ddouble Variance =>
            0.5d * Square(Sigma) * Expm1(1d / Square(Delta)) * (Exp(1d / Square(Delta)) * Cosh(2 * Gamma / Sigma) + 1);

        public override ddouble Skewness {
            get {
                ddouble exp_delta2 = Exp(1d / Square(Delta));

                return Cube(Sigma) * Sqrt(exp_delta2) * Square(Expm1(1d / Square(Delta)))
                        * (exp_delta2 * (exp_delta2 + 2) * Sinh(3 * Gamma / Delta) + 3 * Sinh(2 * Gamma / Delta))
                    / (4 * Cube(Sqrt(Variance)));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble exp_delta2_1 = Exp(1d / Square(Delta));
                ddouble exp_delta2_2 = Square(exp_delta2_1);
                ddouble exp_delta2_3 = Cube(exp_delta2_1);
                ddouble exp_delta2_4 = exp_delta2_2 * exp_delta2_2;

                ddouble k1 = exp_delta2_2 * (exp_delta2_4 + 2 * exp_delta2_3 + 3 * exp_delta2_2 - 3) * Cosh(4 * Gamma / Delta);
                ddouble k2 = 4 * exp_delta2_2 * (exp_delta2_1 + 2) * Cosh(3 * Gamma / Delta);
                ddouble k3 = 3 * (2 * exp_delta2_2 + 1);

                return Pow(Sigma, 4) * Square(Expm1(1d / Square(Delta))) * (k1 + k2 + k3) / (8 * Square(Variance));
            }
        }

        public static JohnsonSUDistribution operator *(JohnsonSUDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma, dist.Mu * k, dist.Sigma * k);
        }

        public static JohnsonSUDistribution operator +(JohnsonSUDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Sigma, dist.Mu + s, dist.Sigma);
        }

        public static JohnsonSUDistribution operator -(JohnsonSUDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Sigma, dist.Mu - s, dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(JohnsonSUDistribution).Name}[gamma={Gamma},delta={Delta},mu={Mu},sigma={Sigma}]";
        }
    }
}
