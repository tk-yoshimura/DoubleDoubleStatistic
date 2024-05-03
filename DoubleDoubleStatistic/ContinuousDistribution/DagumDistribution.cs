using DoubleDouble;
using DoubleDoubleStatistic.Utils;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class DagumDistribution : ContinuousDistribution {

        public ddouble A { get; }
        public ddouble P { get; }

        private readonly ddouble ap, a_inv;

        public DagumDistribution(ddouble a, ddouble p) {
            ValidateShape(a, a => a > 0d);
            ValidateShape(p, p => p > 0d);

            A = a;
            P = p;

            ap = a * p;
            a_inv = 1d / a;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (A == 1d) {
                ddouble pdf = P * Pow(x, P - 1d) / Pow(x + 1d, P + 1d);
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }
            else {
                ddouble xa = Pow(x, A);

                if (xa <= 0d) {
                    return (A < 1d) ? PositiveInfinity : 0d;
                }

                if (IsPositiveInfinity(xa)) {
                    return 0d;
                }

                ddouble pdf = ap * Pow(xa, P) / (x * Pow(xa + 1d, P + 1d));
                pdf = IsFinite(pdf) ? pdf : 0d;

                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble xa = Pow(x, -A), xcp1 = 1d + xa;

            if (interval == Interval.Lower) {
                if (IsNegative(x) || xa <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 1d;
                }

                ddouble cdf = Min(1d, Pow(xcp1, -P));

                return cdf;
            }
            else {
                if (IsNegative(x) || xa <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(x) || IsPositiveInfinity(xcp1)) {
                    return 0d;
                }

                ddouble cdf = Max(0d, 1d - Pow(xcp1, -P));

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

                ddouble x = Pow(Pow(1d / p, 1d / P) - 1d, -a_inv);

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Pow(Pow(1d / (1d - p), 1d / P) - 1d, -a_inv);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            (A > 1d)
            ? P * Beta(P + a_inv, 1d - a_inv)
            : NaN;


        public override ddouble Median =>
            Pow(Pow2(1d / P) - 1d, -a_inv);

        public override ddouble Mode =>
            (P * A > 1d)
            ? Pow((A + 1d) / (P * A - 1d), -a_inv)
            : 0d;

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            (A > 2d)
            ? IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            (A > 3d)
            ? IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            (A > 4d)
            ? IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(DagumDistribution).Name}[a={A},p={P}]";
        }

        public override string Formula => "p(x; a, p) := a * p * x^(a * p - 1) / (1 + x^a)^(p + 1)";
    }
}
