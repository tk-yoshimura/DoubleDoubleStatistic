using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class InverseGaussDistribution : ScalableDistribution<InverseGaussDistribution>,
        IMultiplyOperators<InverseGaussDistribution, ddouble, InverseGaussDistribution> {

        public ddouble Mu { get; }
        public ddouble Lambda { get; }

        private readonly ddouble r;

        public InverseGaussDistribution() : this(mu: 1d, lambda: 1d) { }

        public InverseGaussDistribution(ddouble mu, ddouble lambda) {
            ValidateShape(mu, mu => mu > 0);
            ValidateScale(lambda);

            Mu = mu;
            Lambda = lambda;

            r = Exp(2 * Lambda / Mu);
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

                ddouble cdf = (Erfc(un) + r * Erfc(up)) / 2;
                cdf = Clamp(cdf, 0d, 1d);

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                ddouble cdf = (Erfc(-un) - r * Erfc(up)) / 2;
                cdf = Clamp(cdf, 0d, 1d);

                return cdf;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Mu;

        public override ddouble Median => throw new NotImplementedException();

        public override ddouble Mode {
            get {
                ddouble c = (3d * Mu) / (2d * Lambda);
                return Mu * (Sqrt(1 + Square(c)) - c);
            }
        }

        public override ddouble Variance => Cube(Mu) / Lambda;

        public override ddouble Skewness => 3d * Sqrt(Mu / Lambda);

        public override ddouble Kurtosis => 15d * Mu / Lambda;

        public override ddouble Entropy => throw new NotImplementedException();

        public static InverseGaussDistribution operator *(InverseGaussDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Lambda * k);
        }

        public override string ToString() {
            return $"{typeof(InverseGaussDistribution).Name}[mu={Mu},lambda={Lambda}]";
        }

    }
}
