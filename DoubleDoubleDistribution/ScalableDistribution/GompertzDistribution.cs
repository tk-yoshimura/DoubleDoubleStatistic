using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class GompertzDistribution : ScalableDistribution<GompertzDistribution>,
        IMultiplyOperators<GompertzDistribution, ddouble, GompertzDistribution> {

        public ddouble Eta { get; }
        public ddouble Theta { get; }

        private readonly ddouble theta_inv;

        public GompertzDistribution(ddouble eta, ddouble theta) {
            ValidateShape(eta, eta => eta > 0);
            ValidateScale(theta);

            Eta = eta;
            Theta = theta;

            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble u = x * theta_inv;

            if (IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble pdf = Eta * Exp(Eta * (1d - Exp(u)) + u) * theta_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (interval == Interval.Lower) {
                ddouble cdf = -Expm1(Eta * (1d - Exp(u)));

                return cdf;
            }
            else {
                ddouble cdf = Exp(Eta * (1d - Exp(u)));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Log1p(-Log1p(-p) / Eta) * Theta;

                return x;
            }
            else {
                ddouble x = Log1p(-Log(p) / Eta) * Theta;

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Theta * Exp(Eta) * -Ei(-Eta);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => (Eta >= 1d) ? 0d : -Log(Eta) * Theta;

        public override ddouble Variance => throw new NotImplementedException();

        public override ddouble Skewness => throw new NotImplementedException();

        public override ddouble Kurtosis => throw new NotImplementedException();

        public override ddouble Entropy => throw new NotImplementedException();

        public static GompertzDistribution operator *(GompertzDistribution dist, ddouble k) {
            return new(dist.Eta, dist.Theta * k);
        }

        public override string ToString() {
            return $"{typeof(GompertzDistribution).Name}[k={Eta},theta={Theta}]";
        }
    }
}
