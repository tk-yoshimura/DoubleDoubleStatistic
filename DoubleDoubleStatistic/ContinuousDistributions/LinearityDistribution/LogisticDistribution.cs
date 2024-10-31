using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LogisticDistribution : LinearityDistribution<LogisticDistribution>,
        IAdditionOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        ISubtractionOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        IMultiplyOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        IDivisionOperators<LogisticDistribution, ddouble, LogisticDistribution>,
        IFittableContinuousDistribution<LogisticDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public LogisticDistribution() : this(mu: 0d, sigma: 1d) { }

        public LogisticDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public LogisticDistribution(ddouble mu, ddouble sigma) {
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Mu = mu;
            Sigma = sigma;

            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv, v = Exp(-u);

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble pdf = IsFinite(v) ? (v / Square(1d + v) * sigma_inv) : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble v = Exp(-u), vp1 = v + 1d;

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

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Log(u / (1d - u)) * (double)Sigma + (double)Mu;

            return v;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance =>
            Square(Sigma * Pi) / 3d;

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

        public static (LogisticDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (LogisticDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileLinearFitter<LogisticDistribution>.Fit(new LogisticDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(LogisticDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-u) / (1 + exp(-u))^2 / sigma, u = (x - mu) / sigma";
    }
}
