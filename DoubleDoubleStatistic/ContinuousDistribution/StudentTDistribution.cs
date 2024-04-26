using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class StudentTDistribution : ContinuousDistribution {

        public ddouble Nu { get; }

        private readonly ddouble pdf_norm, nu_inv, nu_half, power;
        private readonly bool is_integer_nu;
        private readonly int n;
        private readonly double zero_thr;

        public bool EnableCDFErrorException { get; set; } = false;

        public StudentTDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            ddouble c = Sqrt(nu * PI);

            pdf_norm = nu < 70d
                ? Gamma((nu + 1) / 2) / (Gamma(nu / 2) * c)
                : Exp(LogGamma((nu + 1) / 2) - LogGamma(nu / 2)) / c;
            nu_inv = 1d / nu;
            nu_half = nu / 2;
            power = -(nu + 1) / 2;
            is_integer_nu = nu <= 1024 && IsInteger(nu);
            n = is_integer_nu ? (int)nu : 0;

            const int zero_thr_log = 710;
            zero_thr = nu < 0.5
                ? double.PositiveInfinity
                : double.Exp((double)(((nu + 1) * Log(nu) + 2 * zero_thr_log) / (2 * nu + 2)));
        }

        public override ddouble PDF(ddouble x) {
            if (Abs(x) >= zero_thr) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = 1 + x * x * nu_inv;
            ddouble v = is_integer_nu ? Pow(Sqrt(u), -(n + 1)) : Pow(u, power);
            ddouble pdf = pdf_norm * v;

            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Upper) {
                return CDF(-x, Interval.Lower);
            }

            ddouble v = Sqrt(x * x + Nu);
            ddouble t = (x + v) / (2 * v);

            ddouble cdf = IncompleteBetaRegularized(t, nu_half, nu_half);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (Abs(p - 0.5d) < 1e-31d) {
                return 0d;
            }

            ddouble x;

            if (p == 0) {
                x = NegativeInfinity;
            }
            else if (p == 1) {
                x = PositiveInfinity;
            }
            else if (Nu == 1) {
                x = Tan(PI * (p - 0.5d));
            }
            else if (Nu == 2) {
                x = 2 * (p - 0.5d) * Sqrt(1d / (2 * p * (1d - p)));
            }
            else if (Nu == 4) {
                ddouble a = 4 * p * (1d - p);
                ddouble q = Cos(Acos(Sqrt(a)) / 3d) / Sqrt(a);

                x = Sign(p - 0.5d) * 2 * Sqrt(q - 1d);
            }
            else {
                ddouble t = InverseIncompleteBeta(p, nu_half, nu_half);
                ddouble u = Sqrt(Nu / (t * (1 - t)));

                x = u * (2 * t - 1) / 2;
            }

            return interval == Interval.Lower ? x : -x;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => Nu > 1d ? 0 : NaN;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Variance => (Nu > 2d)
            ? Nu / (Nu - 2d)
            : (Nu > 1d)
            ? PositiveInfinity
            : NaN;

        public override ddouble Skewness => (Nu > 3d) ? 0d : NaN;

        public override ddouble Kurtosis => (Nu > 4d)
            ? 6d / (Nu - 4d)
            : (Nu > 2d)
            ? PositiveInfinity
            : NaN;

        public override ddouble Entropy =>
            (Nu + 1d) * 0.5d * (Digamma((Nu + 1d) * 0.5d) - Digamma(Nu * 0.5d)) + Log(Sqrt(Nu) * Beta(Nu * 0.5d, 0.5d));

        public override string ToString() {
            return $"{typeof(StudentTDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := (1 + x^2 / nu)^(- (nu + 1) / 2) * gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(pi * nu))";
    }
}
