using DoubleDouble;
using DoubleDoubleIntegrate;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class StudentTDistribution : Distribution {

        public ddouble Nu { get; }

        private readonly ddouble norm, nu_inv, nu_half, power, cdf_norm;
        private readonly bool is_integer_nu;
        private readonly int n;
        private readonly double zero_thr;

        public bool EnableCDFErrorException { get; set; } = false;

        public StudentTDistribution(ddouble nu) {
            ValidateShape(nu);

            this.Nu = nu;

            ddouble c = Sqrt(nu * PI);

            this.norm = nu < 70d
                ? Gamma((nu + 1) / 2) / (Gamma(nu / 2) * c)
                : Exp(LogGamma((nu + 1) / 2) - LogGamma(nu / 2)) / c;
            this.nu_inv = 1d / nu;
            this.nu_half = nu / 2;
            this.power = -(nu + 1) / 2;
            this.cdf_norm = 1d / Beta(nu_half, nu_half);
            this.is_integer_nu = nu <= 1024 && IsInteger(nu);
            this.n = is_integer_nu ? (int)nu : 0;

            const int zero_thr_log = 710;
            this.zero_thr = nu < 0.5
                ? double.PositiveInfinity
                : double.Exp((double)(((nu + 1) * Log(nu) + 2 * zero_thr_log) / (2 * nu + 2)));
        }

        public override ddouble PDF(ddouble x) {
            if (Abs(x) >= zero_thr) {
                return Zero;
            }

            ddouble u = 1 + x * x * nu_inv;
            ddouble v = is_integer_nu ? Pow(Sqrt(u), -(n + 1)) : Pow(u, power);
            ddouble pdf = norm * v;

            pdf = IsFinite(pdf) ? pdf : Zero;

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            if (Abs(x) < 1e-32) {
                return Point5;
            }
            if (x <= -zero_thr) {
                return Zero;
            }
            if (x >= zero_thr) {
                return One;
            }

            if (nu_half <= 64) {
                ddouble u = Sqrt(x * x + Nu), v = (x + u) / Ldexp(u, 1);
                ddouble cdf = IncompleteBeta(v, nu_half, nu_half) * cdf_norm;

                return cdf;
            }
            else {
                return CDFIntegrate(x);
            }
        }

        private ddouble CDFIntegrate(ddouble x) {
            ddouble eps = 1e-27 * norm;

            ddouble u = 1d / (Abs(x) + 1d);

            ddouble error, cdf;

            if (u < 0.5d) {
                ddouble f(ddouble t) {
                    if (IsZero(t)) {
                        return Zero;
                    }

                    ddouble t_inv = 1d / t;
                    ddouble x = (1d - t) * t_inv;

                    ddouble u = 1 + x * x * nu_inv;
                    ddouble pdf = is_integer_nu ? Pow(Sqrt(u), -(n + 1)) : Pow(u, power);

                    ddouble y = pdf * t_inv * t_inv;

                    return y;
                };

                (ddouble value, error) = GaussKronrodIntegral.AdaptiveIntegrate(f, Zero, u, eps, depth: 12);
                value = Max(0d, value) * norm;

                cdf = x < 0d ? value : 1d - value;
            }
            else {
                ddouble f(ddouble t) {
                    ddouble t_inv = 1d / t;
                    ddouble x = (1d - t) * t_inv;

                    ddouble u = 1 + x * x * nu_inv;
                    ddouble pdf = is_integer_nu ? Pow(Sqrt(u), -(n + 1)) : Pow(u, power);

                    ddouble y = pdf * t_inv * t_inv;

                    return y;
                };

                (ddouble value, error) = GaussKronrodIntegral.AdaptiveIntegrate(f, u, One, eps, depth: 12);
                value = Max(0d, value) * norm;

                cdf = x < 0d ? 0.5d - value : 0.5d + value;
            }

            if (EnableCDFErrorException && !(error < eps)) {
                throw new ArithmeticException("CDF integrate not convergence.");
            }

            return cdf;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => (Nu > 1) ? 0 : NaN;
        public override ddouble Median => 0;
        public override ddouble Mode => 0;
        public override ddouble Variance => (Nu > 2) ? (Nu / (Nu - 2)) : ((Nu > 1) ? PositiveInfinity : NaN);
        public override ddouble Skewness => (Nu > 3) ? 0 : NaN;
        public override ddouble Kurtosis => (Nu > 4) ? (6 / (Nu - 4)) : NaN;

        public override ddouble Entropy =>
            (Nu + 1) / 2 * (Digamma((Nu + 1) / 2) - Digamma(Nu / 2)) + Log(Sqrt(Nu) * Beta(Nu / 2, Point5));

        public override string ToString() {
            return $"{typeof(StudentTDistribution).Name}[nu={Nu}]";
        }
    }
}
