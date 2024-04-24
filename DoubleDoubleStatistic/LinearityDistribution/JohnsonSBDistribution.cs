using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class JohnsonSBDistribution : LinearityDistribution<JohnsonSBDistribution>,
        IMultiplyOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        IAdditionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        ISubtractionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution> {

        public ddouble Gamma { get; }
        public ddouble Delta { get; }
        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_inv;

        public JohnsonSBDistribution(ddouble gamma, ddouble delta, ddouble mu, ddouble sigma) {
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

            if (u <= 0d || u >= 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm / (u * (1d - u)) * Exp(-0.5d * Square(Gamma + Delta * Log(u / (1d - u))));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;

            if (u <= 0d) {
                return (interval == Interval.Lower) ? 0d : 1d;
            }
            if (u >= 1d) {
                return (interval == Interval.Lower) ? 1d : 0d;
            }

            ddouble v = Gamma + Delta * Log(u / (1d - u));

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

            ddouble u = 1d / (Exp((Gamma - v) / Delta) + 1);
            ddouble x = Mu + Sigma * u;

            return x;
        }

        public override (ddouble min, ddouble max) Support => (Mu, Mu + Sigma);

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048);

        public override ddouble Mode => throw new NotImplementedException();

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

        public static JohnsonSBDistribution operator +(JohnsonSBDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Sigma, dist.Mu + s, dist.Sigma);
        }

        public static JohnsonSBDistribution operator -(JohnsonSBDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Sigma, dist.Mu - s, dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(JohnsonSBDistribution).Name}[gamma={Gamma},delta={Delta},mu={Mu},sigma={Sigma}]";
        }
    }
}
