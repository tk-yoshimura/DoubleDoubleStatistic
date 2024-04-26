using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class LogNormalDistribution : ContinuousDistribution,
        IMultiplyOperators<LogNormalDistribution, LogNormalDistribution, LogNormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        public LogNormalDistribution() : this(mu: 0, sigma: 1) { }

        public LogNormalDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;

            sigma_sq = sigma * sigma;
            pdf_norm = 1d / (sigma * Sqrt(2 * PI));
            exp_scale = -LbE / (2 * sigma_sq);
            erf_scale = -1d / (Sqrt2 * sigma);
        }

        public override ddouble PDF(ddouble x) {
            if (x <= 0d) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble s = Square(Log(x) - Mu) * exp_scale;

            ddouble pdf = Pow2(s + Log2(pdf_norm / x));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = Ldexp(Erfc((Log(x) - Mu) * erf_scale), -1);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = Ldexp(Erfc((Mu - Log(x)) * erf_scale), -1);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Exp(Mu - Sigma * Sqrt2 * InverseErfc(Ldexp(p, 1)));

                return x;
            }
            else {
                ddouble x = Exp(Mu + Sigma * Sqrt2 * InverseErfc(Ldexp(p, 1)));

                return x;
            }
        }

        public override bool MultiplyClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Exp(Mu + Ldexp(sigma_sq, -1));

        public override ddouble Median => Exp(Mu);

        public override ddouble Mode => Exp(Mu - sigma_sq);

        public override ddouble Variance =>
            Exp(2 * Mu + sigma_sq) * (Exp(sigma_sq) - 1);

        public override ddouble Skewness =>
            Sqrt(Exp(sigma_sq) - 1) * (Exp(sigma_sq) + 2);

        public override ddouble Kurtosis =>
            Exp(4 * sigma_sq) + 2 * Exp(3 * sigma_sq) + 3 * Exp(2 * sigma_sq) - 6d;

        public override ddouble Entropy =>
            (1d + Log(2 * PI * sigma_sq)) / 2 + Mu;

        public static LogNormalDistribution operator *(LogNormalDistribution dist1, LogNormalDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public override string ToString() {
            return $"{typeof(LogNormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-(log(x) - mu)^2 / (2 * sigma^2)) / (x * sigma * sqrt(2 * pi))";
    }
}
