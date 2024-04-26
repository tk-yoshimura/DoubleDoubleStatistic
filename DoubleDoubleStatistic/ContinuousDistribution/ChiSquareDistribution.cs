using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class ChiSquareDistribution : ContinuousDistribution,
        IAdditionOperators<ChiSquareDistribution, ChiSquareDistribution, ChiSquareDistribution> {

        public ddouble Nu { get; }

        private readonly ddouble c, pdf_lognorm;

        public ChiSquareDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            c = Nu * 0.5d - 1d;
            pdf_lognorm = nu * 0.5d + LogGamma(nu * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || IsPositiveInfinity(x)) {
                return 0d;
            }
            if (IsZero(x)) {
                return Nu <= 1 ? PositiveInfinity : Nu <= 2 ? 0.5d : 0d;
            }

            ddouble pdf = Pow2(c * Log2(x) - Ldexp(x, -1) * LbE - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, Ldexp(x, -1));

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, Ldexp(x, -1));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseLowerIncompleteGamma(Nu * 0.5d, p) * 2;

                return x;
            }
            else {
                ddouble x = InverseUpperIncompleteGamma(Nu * 0.5d, p) * 2;

                return x;
            }
        }

        public override bool AdditiveClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Nu;

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => Max(Nu - 2d, 0);

        public override ddouble Variance => Ldexp(Nu, 1);

        public override ddouble Skewness => Sqrt(8d / Nu);

        public override ddouble Kurtosis => 12d / Nu;

        public override ddouble Entropy {
            get {
                ddouble nu_half = Ldexp(Nu, -1);

                return nu_half + Log(2 * Gamma(nu_half)) + (1 - nu_half) * Digamma(nu_half);
            }
        }

        public static ChiSquareDistribution operator +(ChiSquareDistribution dist1, ChiSquareDistribution dist2) {
            return new(dist1.Nu + dist2.Nu);
        }

        public override string ToString() {
            return $"{typeof(ChiSquareDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := x^(nu / 2 - 1) * exp(-x / 2) / (2^(nu / 2) * gamma(nu / 2))";
    }
}
