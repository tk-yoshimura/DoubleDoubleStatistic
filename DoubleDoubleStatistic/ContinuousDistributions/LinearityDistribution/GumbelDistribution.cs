using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class GumbelDistribution : LinearityDistribution<GumbelDistribution>,
        IAdditionOperators<GumbelDistribution, ddouble, GumbelDistribution>,
        ISubtractionOperators<GumbelDistribution, ddouble, GumbelDistribution>,
        IMultiplyOperators<GumbelDistribution, ddouble, GumbelDistribution>,
        IDivisionOperators<GumbelDistribution, ddouble, GumbelDistribution>,
        IFittableContinuousDistribution<GumbelDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public GumbelDistribution() : this(mu: 0d, sigma: 1d) { }

        public GumbelDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public GumbelDistribution(ddouble mu, ddouble sigma) {
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Mu = mu;
            Sigma = sigma;
            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (Mu - x) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsInfinity(u)) {
                return 0d;
            }

            ddouble pdf = Exp(-Exp(u) + u) * sigma_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (Mu - x) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble cdf = Exp(-Exp(u));

                return cdf;
            }
            else {
                ddouble cdf = -Expm1(-Exp(u));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble u = Log(-Log(p));
                ddouble x = Mu - Sigma * u;

                x = IsNaN(x) ? PositiveInfinity : x;

                return x;
            }
            else {
                ddouble u = Log(-Log1p(-p));
                ddouble x = Mu - Sigma * u;

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Log(-double.Log(1d - u));
            double w = (double)Mu - v * (double)Sigma;

            return w;
        }

        public override ddouble Mean =>
            Mu + EulerGamma * Sigma;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Mu;

        public override ddouble Variance =>
            Square(PI * Sigma) / 6d;

        public override ddouble Skewness =>
            12d * Sqrt(6d) * Zeta3 / Cube(PI);

        public override ddouble Kurtosis => (ddouble)12 / 5;

        public override ddouble Entropy => Log(Sigma) + EulerGamma + 1d;

        public static GumbelDistribution operator +(GumbelDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Sigma);
        }

        public static GumbelDistribution operator -(GumbelDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Sigma);
        }

        public static GumbelDistribution operator *(GumbelDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Sigma * k);
        }

        public static GumbelDistribution operator /(GumbelDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Sigma / k);
        }

        public static (GumbelDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (GumbelDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileLinearFitter<GumbelDistribution>.Fit(new GumbelDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(GumbelDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-exp(u) + u) / sigma, u = (x - mu) / sigma";
    }
}
