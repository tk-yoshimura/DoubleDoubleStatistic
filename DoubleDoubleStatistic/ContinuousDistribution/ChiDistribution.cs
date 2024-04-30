using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class ChiDistribution : ContinuousDistribution {

        public ddouble Nu { get; }

        private readonly ddouble c, pdf_lognorm;

        public ChiDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            c = nu * 0.5d - 1d;
            pdf_lognorm = c + LogGamma(nu * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }
            if (IsNegative(x) || IsPositiveInfinity(x)) {
                return 0d;
            }
            if (IsZero(x)) {
                return Nu <= 1 ? Sqrt(2 * RcpPI) : 0d;
            }

            ddouble pdf = Pow2((Nu - 1) * Log2(x) - x * x * LbE * 0.5d - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble x2 = x * x;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, x2 * 0.5d);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, x2 * 0.5d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(InverseLowerIncompleteGamma(Nu * 0.5d, p) * 2d);

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(Nu * 0.5d, p) * 2d);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Sqrt2 * Exp(LogGamma((Nu + 1d) * 0.5d) - LogGamma(Nu * 0.5d));

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Sqrt(Nu - 1d);

        public override ddouble Variance => Nu - Square(Mean);

        public override ddouble Skewness {
            get {
                ddouble variance = Variance;

                return Mean * (1d - 2d * variance) / ExMath.Pow3d2(variance);
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble variance = Variance;

                return 2d * (1d / variance - Mean * Skewness / Sqrt(variance) - 1d);
            }
        }

        public override ddouble Entropy {
            get {
                ddouble nu_half = Nu * 0.5d;

                return LogGamma(nu_half) + (Nu - Ln2 - (Nu - 1d) * Polygamma(0, nu_half)) * 0.5d;
            }
        }

        public override string ToString() {
            return $"{typeof(ChiDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := x^(nu - 1) * exp(-x^2 / 2) / (2^(nu / 2 - 1) * gamma(nu / 2))";
    }
}
