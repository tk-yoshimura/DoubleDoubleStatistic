using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class ExponentialDistribution : ScalableDistribution<ExponentialDistribution>,
        IMultiplyOperators<ExponentialDistribution, ddouble, ExponentialDistribution>,
        IDivisionOperators<ExponentialDistribution, ddouble, ExponentialDistribution> {

        public ddouble Theta { get; }

        private readonly ddouble theta_inv;

        public ExponentialDistribution(ddouble theta) {
            ValidateScale(theta);

            Theta = theta;

            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble pdf = Exp(-u) * theta_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }

                ddouble cdf = -Expm1(-u);

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                ddouble cdf = Exp(-u);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = -Log1p(-p) * Theta;

                return x;
            }
            else {
                ddouble x = -Log(p) * Theta;

                if (IsNegative(x)) {
                    return 0d;
                }

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Theta;

        public override ddouble Median => Ln2 * Theta;

        public override ddouble Mode => 0d;

        public override ddouble Variance => Square(Theta);

        public override ddouble Skewness => 2d;

        public override ddouble Kurtosis => 6d;

        public override ddouble Entropy => 1d + Log(Theta);

        public static ExponentialDistribution operator *(ExponentialDistribution dist, ddouble k) {
            return new(dist.Theta * k);
        }

        public static ExponentialDistribution operator /(ExponentialDistribution dist, ddouble k) {
            return new(dist.Theta / k);
        }

        public override string ToString() {
            return $"{typeof(ExponentialDistribution).Name}[theta={Theta}]";
        }

        public override string Formula => "p(x; theta) := exp(-u) / theta, u = x / theta";
    }
}
