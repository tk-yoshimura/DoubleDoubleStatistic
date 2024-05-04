using DoubleDouble;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    [DebuggerDisplay("{ToString(),nq}")]
    public class JohnsonSBDistribution : LinearityDistribution<JohnsonSBDistribution>,
        IMultiplyOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        IAdditionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        ISubtractionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        IDivisionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution> {

        public ddouble Gamma { get; }
        public ddouble Delta { get; }
        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private readonly ddouble pdf_norm, sigma_inv;

        public JohnsonSBDistribution(ddouble gamma, ddouble delta) : this(gamma, delta, mu: 0d, sigma: 1d) { }

        public JohnsonSBDistribution(ddouble gamma, ddouble delta, ddouble sigma) : this(gamma, delta, mu: 0d, sigma: sigma) { }

        public JohnsonSBDistribution(ddouble gamma, ddouble delta, ddouble mu, ddouble sigma) {
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

            if (u <= 0d || u >= 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm / (u * (1d - u)) * Exp(-0.5d * Square(Gamma + Delta * Log(u / (1d - u))));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (u <= 0d) {
                return (interval == Interval.Lower) ? 0d : 1d;
            }
            if (u >= 1d) {
                return (interval == Interval.Lower) ? 1d : 0d;
            }

            ddouble v = Gamma + Delta * Log(u / (1d - u));

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

            ddouble u = 1d / (Exp((Gamma - v) / Delta) + 1d);
            ddouble x = Mu + Sigma * u;

            return x;
        }

        public override (ddouble min, ddouble max) Support => (Mu, Mu + Sigma);

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble u = 0.5d;

                for (int i = 0; i < 256; i++) {
                    ddouble v = u * (u - 1d);

                    ddouble du = (2d * u - 1d - Delta * (Gamma + Delta * Log(u / (1d - u)))) * v
                        / (2d * v + Delta * Delta);

                    if (!IsFinite(du)) {
                        break;
                    }

                    u = Clamp(u - du, u / 16d, (u + 15d) / 16d);

                    if (Abs(du) <= Abs(u) * 1e-29) {
                        break;
                    }
                }

                ddouble x = mode ??= Mu + u * Sigma;

                return x;
            }
        }

        public override ddouble Median => Quantile(0.5d);

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static JohnsonSBDistribution operator *(JohnsonSBDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma, dist.Mu * k, dist.Sigma * k);
        }

        public static JohnsonSBDistribution operator /(JohnsonSBDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma, dist.Mu / k, dist.Sigma / k);
        }

        public static JohnsonSBDistribution operator +(JohnsonSBDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Sigma, dist.Mu + s, dist.Sigma);
        }

        public static JohnsonSBDistribution operator -(JohnsonSBDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Sigma, dist.Mu - s, dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(JohnsonSBDistribution).Name}[gamma={Gamma},delta={Delta},mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; gamma, delta, mu, sigma) := (delta * exp(-(gamma + delta * log(u / (1 - u)))^2 / 2)) / (sqrt(2 * pi) * u * (1 - u)) / sigma, u = (x - mu) / sigma";
    }
}
