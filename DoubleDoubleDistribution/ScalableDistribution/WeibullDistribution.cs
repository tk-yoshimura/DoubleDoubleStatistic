using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class WeibullDistribution : ScalableDistribution<WeibullDistribution>,
        IMultiplyOperators<WeibullDistribution, ddouble, WeibullDistribution> {

        public ddouble K { get; }
        public ddouble Theta { get; }

        private readonly ddouble k_inv, theta_inv;

        public WeibullDistribution(ddouble k, ddouble theta) {
            ValidateShape(k, k => k > 0d);
            ValidateScale(theta);

            K = k;
            Theta = theta;

            k_inv = 1d / K;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (x == 0d) {
                return (K < 1d) ? PositiveInfinity : (K == 1d) ? theta_inv : 0d;
            }

            ddouble u = x * theta_inv;
            ddouble pdf = Exp(Log(u) * (K - 1d) - Pow(u, K)) * K * theta_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = -Expm1(-Pow(u, K));

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = Exp(-Pow(u, K));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p == 1d) {
                    return PositiveInfinity;
                }

                ddouble x = Pow(-Log1p(-p), k_inv) * Theta;

                return x;
            }
            else {
                if (p == 0d) {
                    return PositiveInfinity;
                }

                ddouble x = Pow(-Log(p), k_inv) * Theta;

                if (IsNegative(x)) {
                    return 0d;
                }

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Theta * Gamma(1d + k_inv);

        public override ddouble Median =>
            Theta * Pow(Ln2, k_inv);

        public override ddouble Mode => (K <= 1d)
            ? 0d
            : Theta * Pow((K - 1d) * k_inv, k_inv);

        public override ddouble Variance =>
            Theta * Theta * (Gamma(1d + 2d * k_inv) - Square(Gamma(1d + k_inv)));

        public override ddouble Skewness {
            get {
                ddouble mu = Mean, var = Variance;

                return (Gamma(1d + 3d * k_inv) * Cube(Theta) - 3d * mu * var - Cube(mu)) / Cube(Sqrt(var));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mu = Mean, var = Variance;

                return (Gamma(1d + 4d * k_inv) * Square(Square(Theta))
                    - 4d * mu * (Gamma(1d + 3d * k_inv) * Cube(Theta) - 3d * mu * var - Cube(mu))
                    - 6d * Square(mu) * var
                    - Square(Square(mu))) /
                    Square(var) - 3d;
            }
        }

        public override ddouble Entropy => 1d + EulerGamma * (1d - k_inv) + Log(Theta * k_inv);

        public static WeibullDistribution operator *(WeibullDistribution dist, ddouble k) {
            return new(dist.K, k * dist.Theta);
        }

        public override string ToString() {
            return $"{typeof(WeibullDistribution).Name}[k={K},theta={Theta}]";
        }
    }
}
