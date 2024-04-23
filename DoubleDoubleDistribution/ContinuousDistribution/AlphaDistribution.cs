using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class AlphaDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm, cdf_norm;

        private static readonly ddouble phi_scale = 1d / Sqrt2;

        public AlphaDistribution() : this(alpha: 1d) { }

        public AlphaDistribution(ddouble alpha) {
            ValidateScale(alpha);

            Alpha = alpha;
            cdf_norm = 1d / Erfc(-alpha * phi_scale);
            pdf_norm = cdf_norm * Sqrt(2d / PI);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = Alpha - 1d / x;

            ddouble pdf = pdf_norm / (x * x) * Exp(-u * u / 2);

            if (IsNaN(pdf)) {
                return 0d;
            }

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = Alpha - 1d / x;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = cdf_norm * Erfc(-u * phi_scale);

                if (IsNaN(cdf)) {
                    return 1d;
                }

                cdf = Min(1d, cdf);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = 1d - cdf_norm * Erfc(-u * phi_scale);

                if (IsNaN(cdf)) {
                    return 0d;
                }

                cdf = Max(0d, cdf);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p, Interval.Lower);
            }

            if (p == 1d) {
                return PositiveInfinity;
            }

            ddouble u = InverseErfc(p / cdf_norm) / phi_scale;
            ddouble x = 1d / Max(0d, Alpha + u);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (Zero, PositiveInfinity);

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => (Sqrt(Alpha * Alpha + 8d) - Alpha) / 4;

        public override ddouble Entropy =>
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(AlphaDistribution).Name}[sigma={Alpha}]";
        }
    }
}
