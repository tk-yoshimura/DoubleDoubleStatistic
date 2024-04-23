using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ChiDistribution : ContinuousDistribution {

        public ddouble Nu { get; }

        private readonly ddouble c, pdf_lognorm;

        public ChiDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            c = Nu * 0.5d - 1d;
            pdf_lognorm = c + LogGamma(nu * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsZero(x)) {
                return Nu <= 1 ? Sqrt(2 * RcpPI) : 0d;
            }

            ddouble pdf = Pow2((Nu - 1) * Log2(x) - Ldexp(x * x, -1) * LbE - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble x2 = x * x;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x2)) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, Ldexp(x * x, -1));

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x2)) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, Ldexp(x * x, -1));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(InverseLowerIncompleteGamma(Nu * 0.5d, p) * 2);

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(Nu * 0.5d, p) * 2);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Sqrt2 * Exp(LogGamma((Nu + 1) * 0.5) - LogGamma(Nu * 0.5));

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => Sqrt(Nu - 1);

        public override ddouble Variance => Nu - Square(Mean);

        public override ddouble Skewness {
            get {
                ddouble variance = Variance;

                return Mean / Cube(Sqrt(variance)) * (1d - 2 * variance);
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
                ddouble k_half = Ldexp(Nu, -1);

                return LogGamma(k_half) + (Nu - Ln2 - (Nu - 1d) * Polygamma(0, k_half)) / 2;
            }
        }

        public override string ToString() {
            return $"{typeof(ChiDistribution).Name}[nu={Nu}]";
        }
    }
}
