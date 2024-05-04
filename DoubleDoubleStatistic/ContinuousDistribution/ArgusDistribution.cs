using DoubleDouble;
using DoubleDoubleStatistic.Utils;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class ArgusDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm, psi, psi_inv, alpha_sq;

        public ArgusDistribution(ddouble alpha) {
            ValidateShape(alpha, alpha => alpha > 0d);

            Alpha = alpha;

            psi = Psi(alpha);
            psi_inv = 1d / psi;
            alpha_sq = alpha * alpha;
            pdf_norm = Cube(alpha) / (Sqrt(2d * PI) * psi);
        }

        private static ddouble Psi(ddouble x) {
            return LowerIncompleteGammaRegularized(1.5d, x * x * 0.5d) * 0.5d;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x < 0d || x > 1d) {
                return 0d;
            }

            ddouble x2 = x * x;

            ddouble pdf = pdf_norm * x * Sqrt(1d - x2) * Exp(-alpha_sq * (1d - x2) * 0.5d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return 1d - CDF(x, Interval.Upper);
            }

            if (x <= 0d) {
                return 1d;
            }
            if (x >= 1d) {
                return 0d;
            }

            ddouble cdf = psi_inv * Psi(Alpha * Sqrt(1d - x * x));

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return Quantile(1d - p, Interval.Upper);
            }

            if (p <= 0d) {
                return 1d;
            }
            if (p >= 1d) {
                return 0d;
            }

            ddouble u = InverseLowerIncompleteGamma(1.5d, psi * p * 2d);
            ddouble v = alpha_sq - 2d * u;

            if (IsNegative(v)) {
                return 0d;
            }

            ddouble x = Sqrt(v) / Alpha;
            x = Min(1d, x);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean =>
            Sqrt(PI / 8d) * Alpha * Exp(-alpha_sq * 0.25d) * BesselI(1, alpha_sq * 0.25d) * psi_inv;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Sqrt((alpha_sq - 2d) + Sqrt(alpha_sq * alpha_sq + 4d)) / (Sqrt2 * Alpha);

        public override ddouble Variance => (1d - 3d / alpha_sq + Alpha * Exp(alpha_sq * -0.5d) / Sqrt(2d * PI) * psi_inv) - Square(Mean);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 8192);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 8192);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(ArgusDistribution).Name}[alpha={Alpha}]";
        }

        public override string Formula => "p(x; alpha) := x * sqrt(1 - x^2) * exp(-alpha^2 * (1 - x^2) / 2) * alpha^3 / (erf(alpha / sqrt(2)) * sqrt(pi / 2) - alpha * exp(-alpha^2 / 2))";
    }
}
