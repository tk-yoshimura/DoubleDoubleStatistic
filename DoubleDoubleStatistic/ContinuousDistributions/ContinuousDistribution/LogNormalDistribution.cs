using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LogNormalDistribution : ContinuousDistribution,
        IMultiplyOperators<LogNormalDistribution, LogNormalDistribution, LogNormalDistribution>,
        IFittableContinuousDistribution<LogNormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        public LogNormalDistribution() : this(mu: 0d, sigma: 1d) { }

        public LogNormalDistribution(ddouble mu, ddouble sigma) {
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Mu = mu;
            Sigma = sigma;

            sigma_sq = sigma * sigma;
            pdf_norm = 1d / (sigma * Sqrt(2d * PI));
            exp_scale = -LbE / (2d * sigma_sq);
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

                ddouble cdf = Erfc((Log(x) - Mu) * erf_scale) * 0.5d;

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = Erfc((Mu - Log(x)) * erf_scale) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Exp(Mu - Sigma * Sqrt2 * InverseErfc(2d * p));

                return x;
            }
            else {
                ddouble x = Exp(Mu + Sigma * Sqrt2 * InverseErfc(2d * p));

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextGaussian();
            double v = double.Exp(u * (double)Sigma + (double)Mu);

            return v;
        }

        public override bool MultiplyClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Exp(Mu + sigma_sq * 0.5d);

        public override ddouble Median => Exp(Mu);

        public override ddouble Mode => Exp(Mu - sigma_sq);

        public override ddouble Variance =>
            Exp(2d * Mu + sigma_sq) * (Exp(sigma_sq) - 1d);

        public override ddouble Skewness =>
            Sqrt(Exp(sigma_sq) - 1d) * (Exp(sigma_sq) + 2d);

        public override ddouble Kurtosis =>
            Exp(4d * sigma_sq) + 2d * Exp(3d * sigma_sq) + 3d * Exp(2d * sigma_sq) - 6d;

        public override ddouble Entropy =>
            (1d + Log(2d * PI * sigma_sq)) / 2d + Mu;

        public static LogNormalDistribution operator *(LogNormalDistribution dist1, LogNormalDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static (LogNormalDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (LogNormalDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble mu = t.u / (1d - Abs(t.u));
                    ddouble sigma = t.v / (1d - t.v);

                    try {
                        LogNormalDistribution dist = new(mu, sigma);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((-0.9999, 1e-10d), (0.9999, 1000d / 1001d)), iter: 64
            );

            try {
                ddouble mu = u / (1d - Abs(u));
                ddouble sigma = v / (1d - v);
                LogNormalDistribution dist = new(mu, sigma);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(LogNormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-(log(x) - mu)^2 / (2 * sigma^2)) / (x * sigma * sqrt(2 * pi))";
    }
}