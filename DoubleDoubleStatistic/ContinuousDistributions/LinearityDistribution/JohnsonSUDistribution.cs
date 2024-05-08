using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class JohnsonSUDistribution : LinearityDistribution<JohnsonSUDistribution>,
        IMultiplyOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution>,
        IDivisionOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution>,
        IAdditionOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution>,
        ISubtractionOperators<JohnsonSUDistribution, ddouble, JohnsonSUDistribution> {

        public ddouble Gamma { get; }
        public ddouble Delta { get; }
        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private readonly ddouble pdf_norm, sigma_inv;

        public JohnsonSUDistribution(ddouble gamma, ddouble delta) : this(gamma, delta, mu: 0d, sigma: 1d) { }

        public JohnsonSUDistribution(ddouble gamma, ddouble delta, ddouble sigma) : this(gamma, delta, mu: 0d, sigma: sigma) { }

        public JohnsonSUDistribution(ddouble gamma, ddouble delta, ddouble mu, ddouble sigma) {
            ValidateShape(gamma, IsFinite);
            ValidateScale(delta);
            ValidateLocation(mu);
            ValidateScale(sigma);

            Gamma = gamma;
            Delta = delta;
            Mu = mu;
            Sigma = sigma;

            pdf_norm = delta / (sigma * Sqrt(2d * PI));
            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble pdf = pdf_norm / Hypot(1d, u) * Exp(-0.5d * Square(Gamma + Delta * Arsinh(u)));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble v = Gamma + Delta * Arsinh(u);

            if (interval == Interval.Lower) {
                ddouble cdf = Erfc(-v * sqrt2_inv) * 0.5d;

                return cdf;
            }
            else {
                ddouble cdf = Erfc(v * sqrt2_inv) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble v = interval == Interval.Lower
                ? -InverseErfc(p * 2d) * Sqrt2
                : InverseErfc(p * 2d) * Sqrt2;

            ddouble u = Sinh((v - Gamma) / Delta);
            ddouble x = Mu + Sigma * u;

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();
            double v = ErrorFunction.InverseErfc(u * 2d) * double.Sqrt(2);
            double w = (double)Mu + (double)Sigma * double.Sinh((v - (double)Gamma) / (double)Delta);

            return w;
        }

        public override ddouble Mean =>
            Mu - Sigma * Exp(1d / (2 * Square(Delta))) * Sinh(Gamma / Delta);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble u = 0d;

                for (int i = 0; i < 256; i++) {
                    ddouble v = u * u + 1d, w = Delta * v * (Gamma + Delta * Arsinh(u)), v_sqrt = Sqrt(v);

                    ddouble du = (w + u * v_sqrt) * v
                        / (2d * u * w + ((Delta * Delta * v) + (v + u * u)) * v_sqrt);

                    if (!IsFinite(du)) {
                        break;
                    }

                    u -= du;

                    if (Abs(du) <= Abs(u) * 1e-29) {
                        break;
                    }
                }

                ddouble x = mode ??= Mu + u * Sigma;

                return x;
            }
        }

        public override ddouble Median =>
            Mu + Sigma * Sinh(-Gamma / Delta);

        public override ddouble Variance =>
            0.5d * Square(Sigma) * Expm1(1d / Square(Delta)) * (Exp(1d / Square(Delta)) * Cosh(2 * Gamma / Delta) + 1);

        public override ddouble Skewness {
            get {
                ddouble omega = Exp(1d / Square(Delta)), omega_m1 = Expm1(1d / Square(Delta));
                ddouble gd = Gamma / Delta;

                return -Sign(gd) * Sqrt(omega * omega_m1 * Square(omega * (2d + omega) * Sinh(3d * gd) + 3 * Sinh(gd))
                    / (2d * Cube(omega * Cosh(2d * gd) + 1d)));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble omega = Exp(1d / Square(Delta));
                ddouble gd = Gamma / Delta;

                return (3d + omega * (6d + omega * ((-3d + omega * omega * (3d + omega * (2d + omega))) * Cosh(4d * gd) + 4d * (2d + omega) * Cosh(2d * gd))))
                    / (2d * Square(omega * Cosh(2d * gd) + 1d)) - 3d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static JohnsonSUDistribution operator *(JohnsonSUDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma, dist.Mu * k, dist.Sigma * k);
        }

        public static JohnsonSUDistribution operator /(JohnsonSUDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma, dist.Mu / k, dist.Sigma / k);
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

        public override string Formula => "p(x; gamma, delta, mu, sigma) := (delta * exp(-(gamma + delta * asinh(u))^2 / 2)) / (sqrt(2 * pi) * sqrt(u^2 + 1)) / sigma, u = (x - mu) / sigma";
    }
}
