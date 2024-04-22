using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class InverseChiDistribution : ContinuousDistribution {

        public ddouble Nu { get; }

        private readonly ddouble c, pdf_lognorm;

        public InverseChiDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            c = Nu * 0.5d + 1d;
            pdf_lognorm = nu * 0.5d - 1d + LogGamma(nu * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble pdf = Pow2(-c * Log2(x) - LbE / (2 * x * x) - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, 1d / (2 * x * x));

                if (IsNaN(cdf)) {
                    return x < Mean ? 0d : 1d;
                }

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, 1d / (2 * x * x));

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

        public override ddouble Mean => Nu > 2d
            ? Exp(LogGamma(Nu - 1d) - LogGamma(Nu * 0.5d)) * Sqrt2
            : NaN;

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => 1d / Sqrt(Nu + 2d);

        public override ddouble Variance => Nu > 4d
            ? 2d / (Nu - 2d) - Square(Mean)
            : NaN;

        public override string ToString() {
            return $"{typeof(InverseChiDistribution).Name}[nu={Nu}]";
        }
    }
}
