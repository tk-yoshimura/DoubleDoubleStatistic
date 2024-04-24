using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class HyperbolicSecantDistribution : ScalableDistribution<HyperbolicSecantDistribution>,
        IMultiplyOperators<HyperbolicSecantDistribution, ddouble, HyperbolicSecantDistribution> {

        public ddouble S { get; }

        private readonly ddouble s_inv, pdf_norm;

        public HyperbolicSecantDistribution() : this(s: 1) { }

        public HyperbolicSecantDistribution(ddouble s) {
            S = s;
            s_inv = 1d / s;
            pdf_norm = s_inv * 0.5d;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * s_inv;

            ddouble pdf = pdf_norm / Cosh(u * PI * 0.5d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble cdf = 2 * RcpPI * Atan(Exp(-Abs(x * s_inv) * PI * 0.5d));
            cdf = Max(cdf, 0d);

            cdf = interval != Interval.Lower ^ IsNegative(x) ? cdf : 1d - cdf;

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble x = Log(TanPI(p * 0.5d)) * RcpPI * 2d;

            x = interval == Interval.Lower ? x : -x;

            x *= S;

            return x;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => 0d;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Variance => S * S;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => 2d;

        public override ddouble Entropy => Log(4d * S);

        public static HyperbolicSecantDistribution operator *(HyperbolicSecantDistribution dist, ddouble k) {
            return new(dist.S * k);
        }

        public override string ToString() {
            return $"{typeof(HyperbolicSecantDistribution).Name}[s={S}]";
        }
    }
}
