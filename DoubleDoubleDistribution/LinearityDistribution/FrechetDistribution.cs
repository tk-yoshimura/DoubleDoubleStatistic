using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class FrechetDistribution : LinearityDistribution<FrechetDistribution>,
        IAdditionOperators<FrechetDistribution, ddouble, FrechetDistribution>,
        ISubtractionOperators<FrechetDistribution, ddouble, FrechetDistribution>,
        IMultiplyOperators<FrechetDistribution, ddouble, FrechetDistribution> {

        public ddouble Alpha { get; }
        public ddouble Mu { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_norm, theta_inv;

        public FrechetDistribution(ddouble alpha, ddouble mu, ddouble theta) {
            ValidateShape(alpha, alpha => alpha > 0);
            ValidateLocation(theta);
            ValidateScale(theta);

            Alpha = alpha;
            Mu = mu;
            Theta = theta;

            pdf_norm = Alpha / Theta;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * theta_inv;
            if (u <= 0d) {
                return 0d;
            }
            if (IsNaN(u)) {
                return NaN;
            }

            ddouble v = Pow(u, Alpha);

            ddouble pdf = this.pdf_norm * Exp(-1d / v) / (u * v);

            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * theta_inv;

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }

                ddouble v = Pow(u, Alpha);

                if (!IsFinite(v)) {
                    return 1d;
                }

                ddouble cdf = Exp(-1d / v);

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }

                ddouble v = Pow(u, Alpha);

                if (!IsFinite(v)) {
                    return 0d;
                }

                ddouble cdf = -Expm1(-1d / v);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return 0d;
                }
                if (p >= 1d) {
                    return PositiveInfinity;
                }

                ddouble v = 1d / -Log(p);
                ddouble u = Pow(v, 1d / Alpha);
                ddouble x = Mu + u * Theta;

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble v = 1d / -Log1p(-p);
                ddouble u = Pow(v, 1d / Alpha);
                ddouble x = Mu + u * Theta;

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (Mu, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d)
            ? Mu + Theta * Gamma(1d - 1d / Alpha)
            : PositiveInfinity;

        public override ddouble Median =>
            Mu + Theta / Pow(Ln2, 1d / Alpha);

        public override ddouble Mode =>
            Mu + Pow(Alpha / (1d + Alpha), 1d / Alpha) * Theta;

        public override ddouble Variance => (Alpha > 2d)
            ? Square(Theta) * (Gamma(1d - 2d / Alpha) - Square(Gamma(1d - 1d / Alpha)))
            : PositiveInfinity;

        public override ddouble Skewness {
            get {
                if (Alpha <= 3d) {
                    return PositiveInfinity;
                }

                ddouble g1 = Gamma(1d - 1d / Alpha), g2 = Gamma(1d - 2d / Alpha), g3 = Gamma(1d - 3d / Alpha);

                return (g3 - 3d * g2 * g1 + 2d * Cube(g1)) / Sqrt(Cube(Variance));
            }
        }

        public override ddouble Kurtosis {
            get {
                if (Alpha <= 4d) {
                    return PositiveInfinity;
                }

                ddouble g1 = Gamma(1d - 1d / Alpha), g2 = Gamma(1d - 2d / Alpha);
                ddouble g3 = Gamma(1d - 3d / Alpha), g4 = Gamma(1d - 4d / Alpha);

                return (g4 - 4d * g3 * g1 + 3d * Square(g2)) / Square(Variance) - 6d;
            }
        }

        public override ddouble Entropy =>
            1d + EulerGamma * (Alpha + 1d) / Alpha + Log(Theta / Alpha);

        public static FrechetDistribution operator +(FrechetDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Theta);
        }

        public static FrechetDistribution operator -(FrechetDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Theta);
        }

        public static FrechetDistribution operator *(FrechetDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Theta * k);
        }

        public override string ToString() {
            return $"{typeof(FrechetDistribution).Name}[alpha={Alpha},mu={Mu},theta={Theta}]";
        }
    }
}
