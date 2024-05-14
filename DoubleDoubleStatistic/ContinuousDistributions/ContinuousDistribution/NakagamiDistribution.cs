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
    public class NakagamiDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<NakagamiDistribution> {

        public ddouble M { get; }
        public ddouble Omega { get; }

        private readonly ddouble pdf_lognorm, momega, omegam;

        private readonly GammaDistribution randam_gen_gamma_dist;

        public NakagamiDistribution(ddouble m, ddouble omega) {
            ValidateShape(m, m => m >= 0.5d);
            ValidateScale(omega);

            M = m;
            Omega = omega;

            momega = m / omega;
            omegam = omega / m;
            pdf_lognorm = -LogGamma(m) * LbE + m * Log2(momega) + 1d;

            randam_gen_gamma_dist = new(kappa: m, theta: omegam);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble pdf = Pow2(Log2(x) * (2d * M - 1d) - momega * x * x * LbE + pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = x * x * momega;

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(M, u);

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(M, u);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(InverseLowerIncompleteGamma(M, p) * omegam);

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(M, p) * omegam);

                return x;
            }
        }

        public override double Sample(Random random) {
            return double.Sqrt(randam_gen_gamma_dist.Sample(random));
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            Exp(LogGamma(M + 0.5d) - LogGamma(M)) * Sqrt(omegam);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode =>
            Sqrt((2d * M - 1d) * Omega / (2d * M));

        public override ddouble Variance =>
            Omega * (1d - Square(Exp(LogGamma(M + 0.5d) - LogGamma(M))) / M);

        public override ddouble Skewness {
            get {
                ddouble mp5 = Exp(LogGamma(M + 0.5d) - LogGamma(M));
                ddouble variance = M - mp5 * mp5;

                return mp5 * (0.5d - 2d * variance) / ExMath.Pow3d2(variance);
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mp5 = Exp(LogGamma(M + 0.5d) - LogGamma(M));
                ddouble variance = M - mp5 * mp5;

                return (M * (4d * M + 1d) - 2d * (2d * M + 1d) * mp5 * mp5) / Square(variance) - 6d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static (NakagamiDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (NakagamiDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = GridMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble m = t.u / (1d - t.u);
                    ddouble omega = t.v / (1d - t.v);

                    try {
                        NakagamiDistribution dist = new(m, omega);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((0.3334d, 1e-10d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble m = u / (1d - u);
                ddouble omega = v / (1d - v);
                NakagamiDistribution dist = new(m, omega);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(NakagamiDistribution).Name}[m={M},omega={Omega}]";
        }

        public override string Formula => "p(x; m, omega) := 2 * x^(2 * m - 1) * exp(- x^2 * m / omega) * (m / omega)^m / gamma(m)";
    }
}
