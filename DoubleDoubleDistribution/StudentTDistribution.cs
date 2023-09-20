using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class StudentTDistribution : Distribution {

        public ddouble Nu { get; }

        private readonly ddouble norm, nu_inv, power;
        private readonly bool is_integer_nu;
        private readonly int n;

        public StudentTDistribution(ddouble nu) {
            ValidateShape(nu);

            this.Nu = nu;

            ddouble c = Sqrt(nu * PI);

            this.norm = nu < 70d
                ? Gamma((nu + 1) / 2) / (Gamma(nu / 2) * c)
                : Exp(LogGamma((nu + 1) / 2) - LogGamma(nu / 2)) / c;
            this.nu_inv = 1d / nu;
            this.power = -(nu + 1) / 2;
            this.is_integer_nu = nu <= 1024 && IsInteger(nu);
            this.n = is_integer_nu ? (int)nu : 0;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = 1 + x * x * nu_inv;
            ddouble v = is_integer_nu ? Pow(Sqrt(u), -(n + 1)) : Pow(u, power);
            ddouble pdf = norm * v;

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            throw new NotImplementedException();
        }

        public override ddouble Quantile(ddouble p) {
            throw new NotImplementedException();
        }

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
