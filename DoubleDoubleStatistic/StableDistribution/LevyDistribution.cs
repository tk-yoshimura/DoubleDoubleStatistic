using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class LevyDistribution : StableDistribution<LevyDistribution>,
        IAdditionOperators<LevyDistribution, LevyDistribution, LevyDistribution>,
        ISubtractionOperators<LevyDistribution, LevyDistribution, LevyDistribution>,
        IAdditionOperators<LevyDistribution, ddouble, LevyDistribution>,
        ISubtractionOperators<LevyDistribution, ddouble, LevyDistribution>,
        IMultiplyOperators<LevyDistribution, ddouble, LevyDistribution>,
        IDivisionOperators<LevyDistribution, ddouble, LevyDistribution> {

        public override ddouble Mu { get; }
        public override ddouble C { get; }

        private readonly ddouble pdf_norm, c_inv;

        public LevyDistribution() : this(mu: 0d, c: 1d) { }

        public LevyDistribution(ddouble c) : this(mu: 0d, c: c) { }

        public LevyDistribution(ddouble mu, ddouble c) {
            ValidateLocation(mu);
            ValidateScale(c);

            Mu = mu;
            C = c;

            c_inv = 1d / c;
            pdf_norm = 1d / (c * Sqrt(2d * PI));
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * c_inv;
            if (IsNaN(u)) {
                return NaN;
            }
            if (u <= 0d) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Exp(-1d / (2d * u)) / ExMath.Pow3d2(u);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }

                ddouble cdf = Erfc(Sqrt(1d / (2d * u)));

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }

                ddouble cdf = Erf(Sqrt(1d / (2d * u)));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Mu + C / (2d * Square(InverseErfc(p)));

                return x;
            }
            else {
                ddouble x = Mu + C / (2d * Square(InverseErf(p)));

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (Mu, PositiveInfinity);

        public override ddouble Mean => PositiveInfinity;

        public override ddouble Median => Mu + C / (2d * Square(InverseErfc(0.5d)));

        public override ddouble Mode => Mu + C / 3d;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy =>
            (1d + 3d * EulerGamma + Log(16d * PI * C * C)) * 0.5d;

        public override ddouble Alpha => 0.5d;

        public override ddouble Beta => 1d;

        public static LevyDistribution operator +(LevyDistribution dist1, LevyDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Square(Sqrt(dist1.C) + Sqrt(dist2.C)));
        }

        public static LevyDistribution operator -(LevyDistribution dist1, LevyDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, Square(Sqrt(dist1.C) + Sqrt(dist2.C)));
        }

        public static LevyDistribution operator +(LevyDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.C);
        }

        public static LevyDistribution operator -(LevyDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.C);
        }

        public static LevyDistribution operator *(LevyDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public static LevyDistribution operator /(LevyDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.C / k);
        }

        public override string ToString() {
            return $"{typeof(LevyDistribution).Name}[mu={Mu},c={C}]";
        }

        public override string Formula => "p(x; mu, c) := exp(-1 / (2 * u)) / u^(3/2) / (sqrt(2 * pi) * c), u = (x - mu) / c";
    }
}
