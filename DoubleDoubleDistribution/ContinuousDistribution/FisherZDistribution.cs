using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class FisherZDistribution : ContinuousDistribution {

        public int D1 { get; }
        public int D2 { get; }

        private readonly ddouble pdf_norm, logd2xd2;

        public FisherZDistribution(int d1, int d2) {
            ValidateShape(d1, d1 => d1 > 0);
            ValidateShape(d2, d2 => d2 > 0);
            _ = checked(d1 + d2);

            D1 = d1;
            D2 = d2;

            pdf_norm = 2 * Exp((Log(D1) * D1 + Log(D2) * D2) * 0.5d - LogBeta(d1 * 0.5d, d2 * 0.5d));
            logd2xd2 = D2 * Log(D2);
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = Exp(D1 * x - Log(D1 * Exp(2 * x) + D2) * (D1 + D2) * 0.5d);
            ddouble pdf = pdf_norm * u;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble d1x = D1 * Exp(2 * x);

            if (interval == Interval.Lower) {
                if (IsNegativeInfinity(x)) {
                    return 0d;
                }
                if (IsPositiveInfinity(x)) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(d1x / (d1x + D2), D1 * 0.5d, D2 * 0.5d);

                return cdf;
            }
            else {
                if (IsNegativeInfinity(x)) {
                    return 1d;
                }
                if (IsPositiveInfinity(x)) {
                    return 0d;
                }

                ddouble cdf = IncompleteBetaRegularized(1d - d1x / (d1x + D2), D2 * 0.5d, D1 * 0.5d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = InverseIncompleteBeta(p, D1 * 0.5d, D2 * 0.5d);
                ddouble x = Log(D2 * u / (D1 * (1d - u))) * 0.5d;

                return x;
            }
            else {
                ddouble u = InverseIncompleteBeta(1d - p, D1 * 0.5d, D2 * 0.5d);
                ddouble x = Log(D2 * u / (D1 * (1d - u))) * 0.5d;

                return x;
            }
        }

        public override ddouble Mean => throw new NotImplementedException();

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => 0d;

        public override ddouble Variance => throw new NotImplementedException();

        public override ddouble Skewness => throw new NotImplementedException();

        public override ddouble Kurtosis => throw new NotImplementedException();

        public override ddouble Entropy => throw new NotImplementedException();

        public override string ToString() {
            return $"{typeof(FisherZDistribution).Name}[d1={D1},d2={D2}]";
        }
    }
}
