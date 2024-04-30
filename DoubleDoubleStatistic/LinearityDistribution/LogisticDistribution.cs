using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class LogisticDistribution : LinearityDistribution<LogisticDistribution>,
        IAdditionOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        ISubtractionOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        IMultiplyOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        IDivisionOperators<LogisticDistribution, ddouble, LogisticDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public LogisticDistribution() : this(mu: 0, sigma: 1) { }

        public LogisticDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;

            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv, v = Exp(-u);

            ddouble pdf = IsFinite(v) ? (v / Square(1d + v) * sigma_inv) : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv, v = Exp(-u), vp1 = v + 1d;

            if (interval == Interval.Lower) {
                ddouble cdf = 1d / vp1;

                return cdf;
            }
            else {
                ddouble cdf = IsFinite(vp1) ? (v / vp1) : 1d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Mu + Sigma * Log(p / (1d - p));

                return x;
            }
            else {
                ddouble x = Mu + Sigma * Log((1d - p) / p);

                return x;
            }
        }

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance =>
            Square(Sigma * PI) / 3d;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => (ddouble)6 / 5;

        public override ddouble Entropy => Log(Sigma) + 2d;

        public static LogisticDistribution operator +(LogisticDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Sigma);
        }

        public static LogisticDistribution operator -(LogisticDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Sigma);
        }

        public static LogisticDistribution operator *(LogisticDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Sigma * k);
        }

        public static LogisticDistribution operator /(LogisticDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(LogisticDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-u) / (1 + exp(-u))^2 / sigma, u = (x - mu) / sigma";
    }
}
