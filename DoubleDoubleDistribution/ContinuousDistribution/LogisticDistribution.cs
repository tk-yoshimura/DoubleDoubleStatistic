using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class LogisticDistribution : ContinuousDistribution,
        IMultiplyOperators<LogisticDistribution, ddouble, LogisticDistribution> {

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
            ddouble u = Exp((Mu - x) * sigma_inv);

            ddouble pdf = u / (Sigma * Square(1d + u));

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = Exp((Mu - x) * sigma_inv);

            if (interval == Interval.Lower) {
                ddouble cdf = 1d / (1d + u);

                return cdf;
            }
            else {
                ddouble cdf = u / (1d + u);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble quantile = Mu + Sigma * Log(p / (1d - p));

                return quantile;
            }
            else {
                ddouble quantile = Mu + Sigma * Log((1d - p) / p);

                return quantile;
            }
        }

        public override bool Scalable => true;
        public override bool Symmetric => true;

        public override ddouble Mean => Mu;
        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;
        public override ddouble Variance => Square(Sigma * PI) / 3d;
        public override ddouble Skewness => 0;
        public override ddouble Kurtosis => (ddouble)6 / 5;

        public override ddouble Entropy => Log(Sigma) + 2;

        public static LogisticDistribution operator *(LogisticDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.Sigma);
        }

        public override string ToString() {
            return $"{typeof(LogisticDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }
    }
}
