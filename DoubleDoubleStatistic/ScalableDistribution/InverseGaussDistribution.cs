using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class InverseGaussDistribution : ScalableDistribution<InverseGaussDistribution>,
        IMultiplyOperators<InverseGaussDistribution, ddouble, InverseGaussDistribution> {

        public ddouble Mu { get; }
        public ddouble Lambda { get; }

        private readonly ddouble r, c;

        private QuantileBuilder quantile_lower_builder = null, quantile_upper_builder = null;

        public InverseGaussDistribution() : this(mu: 1d, lambda: 1d) { }

        public InverseGaussDistribution(ddouble mu, ddouble lambda) {
            ValidateShape(mu, mu => mu > 0);
            ValidateScale(lambda);

            Mu = mu;
            Lambda = lambda;

            c = Lambda / Mu;
            r = Exp(2 * c);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble pdf = Sqrt(Lambda / (2 * PI * Cube(x))) * Exp(-Lambda * Square(x - Mu) / (2 * Square(Mu) * x));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = Sqrt(Lambda / (x * 2)) / Mu;
            ddouble un = u * (Mu - x), up = u * (Mu + x);

            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }
                if (IsPositiveInfinity(x)) {
                    return 1d;
                }

                ddouble cdf = (Erfc(un) + r * Erfc(up)) / 2;
                cdf = Min(cdf, 1d);

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble cdf = (Erfc(-un) - r * Erfc(up)) / 2;
                cdf = Max(cdf, 0d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            ddouble df(ddouble x) {
                ddouble y = c * Exp(-(c * Square(x - 1d)) / (2d * x)) / (Sqrt(2 * c * PI / x) * x * x);
                return y;
            }

            if (interval == Interval.Lower) {
                ddouble f(ddouble x) {
                    ddouble u = Sqrt(c / (x * 2d));
                    ddouble y = (r * Erfc(u * (x + 1d)) + Erfc(u * (1d - x))) / 2;
                    return Min(1d, y);
                }

                this.quantile_lower_builder ??= new QuantileBuilder(0d, 2d, f);

                (ddouble x, ddouble x0, ddouble x1) = quantile_lower_builder.Estimate(p);

                for (int i = 0; i < 8; i++) {
                    ddouble y = f(x), dx = (y - p) / df(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x - dx, x0, x1);

                    if (Abs(dx / x) < 1e-32) {
                        break;
                    }
                }

                x *= Mu;

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                ddouble f(ddouble x) {
                    ddouble u = Sqrt(c / (x * 2d));
                    ddouble y = (-r * Erfc(u * (x + 1d)) + Erfc(u * (x - 1d))) / 2;
                    return Max(0d, y);
                }

                this.quantile_upper_builder ??= new QuantileBuilder(64d, 0d, f);

                (ddouble x, ddouble x0, ddouble x1) = quantile_upper_builder.Estimate(p);

                for (int i = 0; i < 8; i++) {
                    ddouble y = f(x), dx = (y - p) / df(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x + dx, x1, x0);

                    if (Abs(dx / x) < 1e-32) {
                        break;
                    }
                }

                x *= Mu;

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Mu;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode {
            get {
                ddouble c = (3d * Mu) / (2d * Lambda);
                return Mu * (Sqrt(1 + Square(c)) - c);
            }
        }

        public override ddouble Variance => Cube(Mu) / Lambda;

        public override ddouble Skewness => 3d * Sqrt(Mu / Lambda);

        public override ddouble Kurtosis => 15d * Mu / Lambda;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static InverseGaussDistribution operator *(InverseGaussDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Lambda * k);
        }

        public override string ToString() {
            return $"{typeof(InverseGaussDistribution).Name}[mu={Mu},lambda={Lambda}]";
        }

    }
}
