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

            c = nu * 0.5d - 1d;
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
                return Nu <= 1d ? PositiveInfinity : Nu <= 2d ? 0.5d : 0d;
            }

            ddouble pdf = Pow2(c * Log2(x) - x * LbE * 0.5d - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, x * 0.5d);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, x * 0.5d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseLowerIncompleteGamma(Nu * 0.5d, p) * 2d;

                return x;
            }
            else {
                ddouble x = InverseUpperIncompleteGamma(Nu * 0.5d, p) * 2d;

                return x;
            }
        }

        public override bool AdditiveClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Nu;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Max(Nu - 2d, 0d);

        public override ddouble Variance => Nu * 2d;

        public override ddouble Skewness => Sqrt(8d / Nu);

        public override ddouble Kurtosis => 12d / Nu;

        public override ddouble Entropy {
            get {
                ddouble nu_half = Nu * 0.5d;

                return nu_half + Log(2d * Gamma(nu_half)) + (1d - nu_half) * Digamma(nu_half);
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
