using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class InverseGammaDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<InverseGammaDistribution> {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble pdf_lognorm;

        private readonly GammaDistribution randam_gen_gamma_dist;

        public InverseGammaDistribution() : this(alpha: 1d, beta: 1d) { }

        public InverseGammaDistribution(ddouble alpha, ddouble beta) {
            ValidateShape(alpha, alpha => alpha > 0d);
            ValidateShape(beta, beta => beta > 0d);

            Alpha = alpha;
            Beta = beta;

            pdf_lognorm = alpha * Log2(beta) - LogGamma(alpha) * LbE;

            randam_gen_gamma_dist = new(kappa: alpha, theta: 1d / beta);
        }

        public override ddouble PDF(ddouble x) {
            if (x <= 0d) {
                return 0d;
            }

            ddouble pdf = Pow2(-Beta / x * LbE - (Alpha + 1d) * Log2(x) + pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Alpha, Beta / x);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Alpha, Beta / x);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return 0d;
                }

                ddouble x = Beta / InverseUpperIncompleteGamma(Alpha, p);

                return x;
            }
            else {
                if (p >= 1d) {
                    return 0d;
                }

                ddouble x = Beta / InverseLowerIncompleteGamma(Alpha, p);

                return x;
            }
        }

        public override double Sample(Random random) {
            return 1d / randam_gen_gamma_dist.Sample(random);
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Alpha > 1d
            ? Beta / (Alpha - 1d)
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode =>
            Beta / (Alpha + 1d);

        public override ddouble Variance => Alpha > 2d
            ? Square(Beta / (Alpha - 1d)) / (Alpha - 2d)
            : NaN;

        public override ddouble Skewness => Alpha > 3d
            ? 4d * Sqrt(Alpha - 2d) / (Alpha - 3d)
            : NaN;

        public override ddouble Kurtosis => Alpha > 4d
            ? (30d * Alpha - 66d) / ((Alpha - 3d) * (Alpha - 4d))
            : NaN;

        public override ddouble Entropy => Alpha + Log(Beta) + LogGamma(Alpha) - (Alpha + 1d) * Digamma(Alpha);

        public static (InverseGammaDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (InverseGammaDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = GridMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble alpha = t.u / (1d - t.u);
                    ddouble beta = t.v / (1d - t.v);

                    try {
                        InverseGammaDistribution dist = new(alpha, beta);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-10d, 1e-10d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble alpha = u / (1d - u);
                ddouble beta = v / (1d - v);
                InverseGammaDistribution dist = new(alpha, beta);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(InverseGammaDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := x^(- alpha - 1) * exp(-beta / x) * beta^alpha / gamma(alpha)";
    }
}
