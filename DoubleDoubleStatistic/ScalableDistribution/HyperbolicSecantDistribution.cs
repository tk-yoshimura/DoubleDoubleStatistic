using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class HyperbolicSecantDistribution : ScalableDistribution<HyperbolicSecantDistribution>,
        IMultiplyOperators<HyperbolicSecantDistribution, ddouble, HyperbolicSecantDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv, pdf_norm;

        public HyperbolicSecantDistribution() : this(sigma: 1) { }

        public HyperbolicSecantDistribution(ddouble sigma) {
            Sigma = sigma;
            sigma_inv = 1d / sigma;
            pdf_norm = sigma_inv * 0.5d;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            ddouble pdf = pdf_norm / Cosh(u * PI * 0.5d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;

            ddouble cdf = 2 * RcpPI * Atan(Exp(-Abs(u) * PI * 0.5d));
            cdf = Max(cdf, 0d);

            cdf = interval != Interval.Lower ^ IsNegative(x) ? cdf : 1d - cdf;

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble x = Log(TanPI(p * 0.5d)) * RcpPI * Sigma * 2d;

            if (IsNaN(x)) {
                x = p < 0.5d ? NegativeInfinity : PositiveInfinity;
            }

            x = interval == Interval.Lower ? x : -x;

            return x;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => 0d;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Variance => Sigma * Sigma;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => 2d;

        public override ddouble Entropy => Log(4d * Sigma);

        public static HyperbolicSecantDistribution operator *(HyperbolicSecantDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public override string ToString() {
            return $"{typeof(HyperbolicSecantDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := 1 / cosh(u * pi / 2) / (2 * sigma), u = x / sigma";
    }
}
