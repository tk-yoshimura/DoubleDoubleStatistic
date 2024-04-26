using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class NakagamiDistribution : ContinuousDistribution {

        public ddouble M { get; }
        public ddouble Omega { get; }

        private readonly ddouble pdf_lognorm, momega, omegam;

        public NakagamiDistribution(ddouble m, ddouble omega) {
            ValidateShape(m, m => m >= 0.5d);
            ValidateScale(omega);

            M = m;
            Omega = omega;

            momega = m / omega;
            omegam = omega / m;
            pdf_lognorm = -LogGamma(m) * LbE + m * Log2(momega) + 1d;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x)) {
                return 0d;
            }

            ddouble pdf = Pow2(Log2(x) * (2 * M - 1d) - momega * x * x * LbE + pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = LowerIncompleteGammaRegularized(M, x * x * momega);

                return cdf;
            }
            else {
                ddouble cdf = UpperIncompleteGammaRegularized(M, x * x * momega);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(InverseLowerIncompleteGamma(M, p) * omegam);

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(M, p) * omegam);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Exp(LogGamma(M + 0.5d) - LogGamma(M)) * Sqrt(omegam);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode =>
            Sqrt((2 * M - 1d) * Omega / (2 * M));

        public override ddouble Variance =>
            Omega * (1d - Square(Exp(LogGamma(M + 0.5d) - LogGamma(M))) / M);

        public override ddouble Skewness {
            get {
                ddouble mp5 = Exp(LogGamma(M + 0.5d) - LogGamma(M));
                ddouble w = M - mp5 * mp5;

                return mp5 * (0.5d - 2 * w) / Cube(Sqrt(w));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mp5 = Exp(LogGamma(M + 0.5d) - LogGamma(M));

                return (M * (4 * M + 1d) - 2d * (2 * M + 1) * mp5 * mp5) / Square(M - mp5 * mp5) - 6d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public override string ToString() {
            return $"{typeof(NakagamiDistribution).Name}[m={M},omega={Omega}]";
        }

        public override string Formula => "p(x; m, omega) := 2 * x^(2 * m - 1) * exp(- x^2 * m / omega) * (m / omega)^m / gamma(m)";
    }
}
