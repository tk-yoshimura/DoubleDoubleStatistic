using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class NoncentralBetaDistribution : ContinuousDistribution {

        public ddouble Alpha { get; }
        public ddouble Beta { get; }
        public ddouble Lambda { get; }

        private const int series_maxiter = 8192;
        private const int cache_samples = 1024;

        private readonly ddouble logbeta;
        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null;

        private readonly NoncentralChiSquareDistribution randam_gen_noncchisq_dist;
        private readonly ChiSquareDistribution randam_gen_chisq_dist;

        public NoncentralBetaDistribution(ddouble alpha, ddouble beta, ddouble lambda) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateShape(nameof(beta), ParamAssert.IsFinitePositive(beta));
            ParamAssert.ValidateNonCentricity(nameof(lambda), ParamAssert.IsFinitePositive(lambda));

            Alpha = alpha;
            Beta = beta;
            Lambda = lambda;

            logbeta = LogBeta(Alpha, Beta) * LbE;

            randam_gen_noncchisq_dist = new(alpha * 2d, lambda);
            randam_gen_chisq_dist = new(beta * 2d);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d || x >= 1d) {
                return 0d;
            }

            ddouble ds = Pow2(Log2(x) * (Alpha - 1d) + Log2(1d - x) * (Beta - 1d) - logbeta);
            ddouble s = ds;

            ddouble exp = Lambda * LbE * -0.5d;

            for (int i = 1; i <= series_maxiter; i++) {
                ds *= Lambda * x * (Alpha + Beta + i - 1) / ((Alpha + i - 1) * (2 * i));

                s += ds;

                if (ds <= s * 1e-30 || !IsFinite(ds)) {
                    break;
                }
                if (i >= series_maxiter) {
                    throw new ArithmeticException($"{this}: pdf calculation not convergence.");
                }

                // Overflow avoidance by rescaling
                if (double.ILogB((double)s) >= 16 || double.ILogB((double)ds) >= 16) {
                    (s, ds) = (Ldexp(s, -16), Ldexp(ds, -16));
                    exp += 16d;
                }
            }

            ddouble pdf = Pow2(exp) * s;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return 1d - CDF(x, Interval.Lower);
            }

            if (x <= 0d) {
                return 0d;
            }
            if (x >= 1d) {
                return 1d;
            }

            ddouble beta0 = IncompleteBetaRegularized(x, Alpha, Beta);
            ddouble beta1 = IncompleteBetaRegularized(x, Alpha + 1d, Beta);

            ddouble lambda_half = Lambda * 0.5d;
            ddouble c = 1d;
            ddouble s = beta0;

            ddouble ai = Alpha, ci = ai + Beta - 1d;
            ddouble exp = -lambda_half * LbE;

            for (int j = 1; j <= series_maxiter; j++) {
                c *= lambda_half / j;

                ddouble ds = c * beta1;

                s += ds;

                if (ds <= s * 1e-30 || !IsFinite(ds)) {
                    break;
                }
                if (j >= series_maxiter) {
                    throw new ArithmeticException($"{this}: cdf calculation not convergence.");
                }

                ai += 1d;
                ci += 1d;

                if (beta1 > 1e-12) {
                    (beta1, beta0) = (((ai + ci * x) * beta1 - ci * x * beta0) / ai, beta1);
                }
                else {
                    Debug.WriteLine(
                        "reset recurr incomp.beta: \n" +
                        $"{((ai + ci * x) * beta1 - ci * x * beta0) / ai} -> {IncompleteBetaRegularized(x, Alpha + (j + 1), Beta)}"
                    );

                    (beta1, beta0) = (IncompleteBetaRegularized(x, Alpha + (j + 1), Beta), beta1);
                }

                // Overflow avoidance by rescaling
                if (double.ILogB((double)s) >= 16 || double.ILogB((double)c) >= 16) {
                    (s, c) = (Ldexp(s, -16), Ldexp(c, -16));
                    exp += 16d;
                }
            }

            ddouble cdf = Min(1d, Pow2(exp) * s);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            if (interval == Interval.Lower) {
                this.quantile_lower_builder ??= new QuantileBuilder(
                    0d, 1d,
                    t => CDF(t, Interval.Lower),
                    cache_samples
                );

                (ddouble x, ddouble x0, ddouble x1) = quantile_lower_builder.Estimate(p);

                for (int i = 0; i < 64; i++) {
                    ddouble y = CDF(x, Interval.Lower), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    if (x0 > 0d) {
                        x = Clamp(x - dx, x0, x1);
                    }
                    else {
                        x = Clamp(x - dx, x / 16d, x1);
                    }

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                this.quantile_upper_builder ??= new QuantileBuilder(
                    1d, 0d,
                    t => CDF(t, Interval.Upper),
                    cache_samples
                );

                (ddouble x, ddouble x0, ddouble x1) = quantile_upper_builder.Estimate(p);

                for (int i = 0; i < 64; i++) {
                    ddouble y = CDF(x, Interval.Upper), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x + dx, x1, x0);

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }
        }

        public override double Sample(Random random) {
            double t = randam_gen_noncchisq_dist.Sample(random);
            double s = randam_gen_chisq_dist.Sample(random);

            double r = t / double.Max(double.Epsilon, t + s);

            return r;
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048);

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                mode ??= BisectionMinimizeSearch1D.Search(t => -PDF(t), (0d, 1d), iter: 1024);

                return mode.Value;
            }
        }

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(NoncentralBetaDistribution).Name}[alpha={Alpha},beta={Beta},lambda={Lambda}]";
        }

        public override string Formula =>
            "p(x; n, m, lambda) := exp(-lambda / 2) * sum((lambda / 2)^k / k! * x^(alpha - 1 + k) * (1 - x)^(beta - 1) / beta(alpha + k, beta), k, 0, inf)";
    }
}
