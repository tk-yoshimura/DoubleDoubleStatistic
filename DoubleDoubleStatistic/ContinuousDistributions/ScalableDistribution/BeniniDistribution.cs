using DoubleDouble;
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
    public class BeniniDistribution : ScalableDistribution<BeniniDistribution>,
        IMultiplyOperators<BeniniDistribution, ddouble, BeniniDistribution>,
        IDivisionOperators<BeniniDistribution, ddouble, BeniniDistribution>,
        IFittableContinuousDistribution<BeniniDistribution> {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public BeniniDistribution(ddouble beta, ddouble sigma)
            : this(alpha: 0d, beta: beta, sigma: sigma) { }

        public BeniniDistribution(ddouble alpha, ddouble beta, ddouble sigma) {
            ParamAssert.ValidateShape(nameof(alpha), alpha >= 0 && IsFinite(alpha));
            ParamAssert.ValidateShape(nameof(beta), ParamAssert.IsFinitePositive(beta));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Alpha = alpha;
            Beta = beta;
            Sigma = sigma;

            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (u < 1d || IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble lnu = Log(u), u_inv = 1d / u;

            ddouble pdf = Exp(-lnu * (Alpha + lnu * Beta)) * (u_inv * (Alpha + 2d * Beta * lnu)) * sigma_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;
            ddouble lnu = Log(u);

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (u < 1d) {
                    return 0d;
                }

                if (IsPositiveInfinity(u)) {
                    return 1d;
                }

                ddouble cdf = -Expm1(-lnu * (Alpha + lnu * Beta));

                return cdf;
            }
            else {
                if (u < 1d) {
                    return 1d;
                }

                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = Exp(-lnu * (Alpha + lnu * Beta));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p == 0d) {
                return interval == Interval.Lower ? Sigma : PositiveInfinity;
            }
            if (p == 1d) {
                return interval == Interval.Lower ? PositiveInfinity : Sigma;
            }

            ddouble lnu;

            if (Alpha > 0d) {
                lnu = interval == Interval.Lower
                    ? (Sqrt(Max(0d, Alpha * Alpha - 4d * Beta * Log1p(-p))) - Alpha) / (2d * Beta)
                    : (Sqrt(Max(0d, Alpha * Alpha - 4d * Beta * Log(p))) - Alpha) / (2d * Beta);
            }
            else {
                lnu = interval == Interval.Lower
                    ? Sqrt(Max(0d, -Log1p(-p) / Beta))
                    : Sqrt(Max(0d, -Log(p) / Beta));
            }

            ddouble x = Exp(lnu) * Sigma;

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double alpha = (double)Alpha, beta = (double)Beta;

            double y;
            if (Alpha > 0d) {
                y = (double.Sqrt(alpha * alpha - 4d * beta * double.Log(u)) - alpha) / (2d * beta);
            }
            else {
                y = double.Sqrt(-double.Log(u) / beta);
            }

            double r = double.Exp(y) * (double)Sigma;

            return r;
        }

        public override (ddouble min, ddouble max) Support => (Sigma, PositiveInfinity);

        public override ddouble Mean =>
            Sigma *
            (Sqrt(Pi) * Exp(Square(Alpha - 1d) / (4d * Beta)) * Erfc((Alpha - 1d) / (2d * Sqrt(Beta))) / (2d * Sqrt(Beta)) + 1d);

        public override ddouble Mode =>
            Max(1d, Exp((Sqrt(8d * Beta + 1d) - 2d * Alpha - 1d) / (4d * Beta))) * Sigma;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Variance =>
            Sigma * Sigma *
            (
                -Pi / (4 * Beta) * Exp(Square(Alpha - 1d) / (2d * Beta)) * Square(Erfc((Alpha - 1d) / (2d * Sqrt(Beta))))
                - Sqrt(Pi) / Sqrt(Beta) * (
                    Exp(Square(Alpha - 1d) / (4d * Beta)) * Erfc((Alpha - 1d) / (2d * Sqrt(Beta)))
                    - Exp(Square(Alpha - 2d) / (4d * Beta)) * Erfc((Alpha - 2d) / (2d * Sqrt(Beta)))
                )
            );

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(new BeniniDistribution(Alpha, Beta, 1d), eps: 1e-28, discontinue_eval_points: 8192);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(new BeniniDistribution(Alpha, Beta, 1d), eps: 1e-28, discontinue_eval_points: 8192);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(new BeniniDistribution(Alpha, Beta, 1d), eps: 1e-28, discontinue_eval_points: 8192) + Log(Sigma);

        public static (BeniniDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BeniniDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble alpha = t.u / (1d - t.u);
                    ddouble beta = t.v / (1d - t.v);

                    try {
                        BeniniDistribution dist = new(alpha, beta, 1d);
                        return QuantileScaleFitter<BeniniDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((0, 1e-10d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble alpha = u / (1d - u);
                ddouble beta = v / (1d - v);
                BeniniDistribution dist = new(alpha, beta, 1d);
                return QuantileScaleFitter<BeniniDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public static BeniniDistribution operator *(BeniniDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Beta, dist.Sigma * k);
        }

        public static BeniniDistribution operator /(BeniniDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Beta, dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(BeniniDistribution).Name}[alpha={Alpha},beta={Beta},sigma={Sigma}]";
        }

        public override string Formula => "p(x, alpha, beta, sigma) := exp(-(log(u) * alpha + log(u)^2 * beta)) * ((alpha + 2 * beta * log(u)) / u) / sigma, u = x / sigma";
    }
}
