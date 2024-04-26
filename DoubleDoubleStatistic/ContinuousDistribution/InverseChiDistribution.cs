using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class InverseChiDistribution : ContinuousDistribution {

        public ddouble Nu { get; }

        private readonly ddouble c, pdf_lognorm;

        public InverseChiDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            c = Nu + 1d;
            pdf_lognorm = nu * 0.5d - 1d + LogGamma(nu * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble pdf = Pow2(-c * Log2(x) - LbE / (2 * x * x) - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = 1d / (2 * x * x);

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, u);

                if (IsNaN(cdf)) {
                    return x < Mean ? 0d : 1d;
                }

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, u);

                if (IsNaN(cdf)) {
                    return x < Mean ? 1d : 0d;
                }

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(0.5d / InverseUpperIncompleteGamma(Nu * 0.5d, p));

                return x;
            }
            else {
                ddouble x = Sqrt(0.5d / InverseLowerIncompleteGamma(Nu * 0.5d, p));

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => 1d / Sqrt(Nu + 1d);

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            (Nu > 1d)
            ? IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            (Nu > 2d)
            ? IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            (Nu > 3d)
            ? IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            (Nu > 4d)
            ? IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(InverseChiDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := x^(-nu - 1) * exp(-1 / (2 * x^2)) / (2^(nu / 2 - 1) * gamma(nu / 2))";
    }
}
