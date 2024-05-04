using DoubleDouble;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    [DebuggerDisplay("{ToString(),nq}")]
    public class SnedecorFDistribution : ContinuousDistribution {

        public ddouble N { get; }
        public ddouble M { get; }

        private readonly ddouble pdf_lognorm;

        public SnedecorFDistribution(ddouble n, ddouble m) {
            ValidateShape(n, n => n > 0d);
            ValidateShape(m, m => m > 0d);

            N = n;
            M = m;

            pdf_lognorm = m * Log2(m) * 0.5d - LogBeta(n * 0.5d, m * 0.5d) * LbE;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x)) {
                return 0d;
            }

            ddouble u = N * x;

            if (u <= 0d) {
                return (N < 2d) ? PositiveInfinity : (N == 2d) ? 1d : 0d;
            }

            if (IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble pdf = Pow2((N * Log2(u) - (N + M) * Log2(u + M)) * 0.5d + pdf_lognorm) / x;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = N * x, upm = u + M;

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(upm)) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(u / upm, N * 0.5d, M * 0.5d);

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(upm)) {
                    return 0d;
                }

                ddouble cdf = IncompleteBetaRegularized(M / upm, M * 0.5d, N * 0.5d);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = InverseIncompleteBeta(p, N * 0.5d, M * 0.5d);
                ddouble x = M * u / (N * (1d - u));

                return x;
            }
            else {
                ddouble u = InverseIncompleteBeta(p, M * 0.5d, N * 0.5d);
                ddouble x = M * (1d - u) / (N * u);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => (M > 2d)
            ? M / (M - 2d)
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => (N > 2d)
            ? (N - 2d) * M / (N * (M + 2d))
            : NaN;

        public override ddouble Variance => (M > 4d)
            ? 2d * M * M * (N + M - 2d) / (N * (M - 2d) * (M - 2d) * (M - 4d))
            : NaN;

        public override ddouble Skewness => (M > 6d)
            ? (2d * N + M - 2d) * Sqrt(8d * (M - 4d)) / ((M - 6d) * Sqrt(N * (N + M - 2d)))
            : NaN;

        public override ddouble Kurtosis => (M > 8d)
            ? (12d * (N * (5d * M - 22d) * (N + M - 2d) + (M - 4d) * (M - 2d) * (M - 2d))
                / (N * (M - 6d) * (M - 8d) * (N + M - 2d)))
            : NaN;

        public override ddouble Entropy =>
            LogGamma(N * 0.5d) + LogGamma(M * 0.5d) - LogGamma((N + M) * 0.5d)
            + (1d - N * 0.5d) * Digamma(N * 0.5d) - (1d + M * 0.5d) * Digamma(M * 0.5d)
            + (N + M) * 0.5d * Digamma((N + M) * 0.5d) - Log(N / M);

        public override string ToString() {
            return $"{typeof(SnedecorFDistribution).Name}[n={N},m={M}]";
        }

        public override string Formula => "p(x; n, m) := x^(n / 2 - 1) * (m + n * x)^((- m - n) / 2) * m^(m / 2) * n^(n / 2) / beta(n / 2, m / 2)";
    }
}
