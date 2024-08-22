using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BetaDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<BetaDistribution> {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }

        private readonly ddouble pdf_norm;

        private readonly GammaDistribution randam_gen_gamma_dist_alpha;
        private readonly GammaDistribution randam_gen_gamma_dist_beta;

        public BetaDistribution(ddouble alpha, ddouble beta) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateShape(nameof(beta), ParamAssert.IsFinitePositive(beta));

            Alpha = alpha;
            Beta = beta;

            pdf_norm = 1d / Beta(alpha, beta);

            randam_gen_gamma_dist_alpha = new GammaDistribution(kappa: alpha);
            randam_gen_gamma_dist_beta = new GammaDistribution(kappa: beta);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x < 0d || x > 1d) {
                return 0d;
            }

            if (Alpha != 1d && Beta != 1d) {
                ddouble pdf = pdf_norm * Pow(x, Alpha - 1d) * Pow1p(-x, Beta - 1d);

                return pdf;
            }
            else if (Alpha != 1d) {
                ddouble pdf = pdf_norm * Pow(x, Alpha - 1d);

                return pdf;
            }
            else if (Beta != 1d) {
                ddouble pdf = pdf_norm * Pow1p(-x, Beta - 1d);

                return pdf;
            }
            else {
                return 1d;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }
                if (x >= 1d) {
                    return 1d;
                }

                ddouble cdf = IncompleteBetaRegularized(x, Alpha, Beta);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }
                if (x >= 1d) {
                    return 0d;
                }

                ddouble cdf = IncompleteBetaRegularized(1d - x, Beta, Alpha);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseIncompleteBeta(p, Alpha, Beta);

                return x;
            }
            else {
                ddouble x = InverseIncompleteBeta(1d - p, Alpha, Beta);

                return x;
            }
        }

        public override double Sample(Random random) {
            double r1 = randam_gen_gamma_dist_alpha.Sample(random);
            double r2 = randam_gen_gamma_dist_beta.Sample(random);

            return r1 / double.Max(r1 + r2, double.Epsilon);
        }

        public override bool Symmetric => Alpha == Beta;

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean => Alpha / (Alpha + Beta);

        public override ddouble Median => InverseIncompleteBeta(0.5d, Alpha, Beta);

        public override ddouble Mode =>
            Alpha > 1d && Beta > 1d ? (Alpha - 1d) / (Alpha + Beta - 2d) :
            Alpha <= 1d && Beta >= 1d && Alpha != Beta ? 0d :
            Alpha >= 1d && Beta <= 1d && Alpha != Beta ? 1d :
            NaN;

        public override ddouble Variance =>
            Alpha * Beta / (Square(Alpha + Beta) * (Alpha + Beta + 1d));

        public override ddouble Skewness =>
            2d * (Beta - Alpha) * Sqrt(Alpha + Beta + 1d) / ((Alpha + Beta + 2d) * Sqrt(Alpha * Beta));

        public override ddouble Kurtosis =>
            6d * (Square(Alpha - Beta) * (Alpha + Beta + 1d) - Alpha * Beta * (Alpha + Beta + 2d)) / (Alpha * Beta * (Alpha + Beta + 2d) * (Alpha + Beta + 3d));

        public override ddouble Entropy =>
            LogBeta(Alpha, Beta) - (Alpha - 1d) * Digamma(Alpha) - (Beta - 1d) * Digamma(Beta) + (Alpha + Beta - 2d) * Digamma(Alpha + Beta);

        public static (BetaDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (BetaDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble alpha = t.u / (1d - t.u);
                    ddouble beta = t.v / (1d - t.v);

                    try {
                        BetaDistribution dist = new(alpha, beta);
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
                BetaDistribution dist = new(alpha, beta);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(BetaDistribution).Name}[alpha={Alpha},beta={Beta}]";
        }

        public override string Formula => "p(x; alpha, beta) := x^(alpha - 1) * (1 - x)^(beta - 1) / beta(alpha, beta)";
    }
}
