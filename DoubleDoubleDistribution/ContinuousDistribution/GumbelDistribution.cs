using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class ExtremeValueDistribution : ContinuousDistribution,
        IAdditionOperators<ExtremeValueDistribution, ddouble, ExtremeValueDistribution>,
        IMultiplyOperators<ExtremeValueDistribution, ddouble, ExtremeValueDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        public ExtremeValueDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateShape(sigma, beta => beta > 0);

            Mu = mu;
            Sigma = sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (Mu - x) / Sigma;

            ddouble pdf = Exp(-Exp(u) + u) / Sigma;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (Mu - x) / Sigma;

            if (interval == Interval.Lower) {
                ddouble cdf = Exp(-Exp(u));

                return cdf;
            }
            else {
                ddouble cdf = -Expm1(-Exp(u));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = Log(-Log(p));
                ddouble quantile = Mu - Sigma * u;

                quantile = IsNaN(quantile) ? PositiveInfinity : quantile;

                return quantile;
            }
            else {
                ddouble u = Log(-Log1p(-p));
                ddouble quantile = Mu - Sigma * u;

                return quantile;
            }
        }

        public override bool Scalable => true;
        public override bool Shiftable => true;

        public override ddouble Mean =>
            Mu + EulerGamma * Sigma;

        public override ddouble Median => Quantile(0.5);

        public override ddouble Mode => Mu;

        public override ddouble Variance =>
            Square(PI * Sigma) / 6d;

        public override ddouble Skewness =>
            12d * Sqrt(6d) * Zeta3 / Cube(PI);

        public override ddouble Kurtosis => (ddouble)12 / 5;

        public override ddouble Entropy => Log(Sigma) + EulerGamma + 1d;

        public static ExtremeValueDistribution operator +(ExtremeValueDistribution dist, ddouble s) {
            return new(s + dist.Mu, dist.Sigma);
        }

        public static ExtremeValueDistribution operator *(ExtremeValueDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(ExtremeValueDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }
    }
}
