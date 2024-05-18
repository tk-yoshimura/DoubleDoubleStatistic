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
    public class GammaDistribution : ScalableDistribution<GammaDistribution>,
        IMultiplyOperators<GammaDistribution, ddouble, GammaDistribution>,
        IDivisionOperators<GammaDistribution, ddouble, GammaDistribution>,
        IFittableContinuousDistribution<GammaDistribution> {

        public ddouble Kappa { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_lognorm, theta_inv;

        private bool randam_gen_param_initialized = false;
        private (double c1, double c2, double c3) randam_gen_param;

        public GammaDistribution(ddouble kappa) : this(kappa, theta: 1d) { }

        public GammaDistribution(ddouble kappa, ddouble theta) {
            ParamAssert.ValidateShape(nameof(kappa), ParamAssert.IsFinitePositive(kappa));
            ParamAssert.ValidateScale(nameof(theta), ParamAssert.IsFinitePositive(theta));

            Kappa = kappa;
            Theta = theta;

            pdf_lognorm = LogGamma(kappa) * LbE;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u) || IsPositiveInfinity(u)) {
                return 0d;
            }

            if (IsZero(u)) {
                return Kappa < 1d ? PositiveInfinity : Kappa == 1d ? theta_inv : 0d;
            }

            ddouble pdf = Pow2((Kappa - 1d) * Log2(u) - u * LbE - pdf_lognorm) * theta_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Kappa, u);

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Kappa, u);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = InverseLowerIncompleteGamma(Kappa, p) * Theta;

                return x;
            }
            else {
                ddouble x = InverseUpperIncompleteGamma(Kappa, p) * Theta;

                return x;
            }
        }

        public override double Sample(Random random) {
            if (!randam_gen_param_initialized) {
                SetupRandamGenParams();
                randam_gen_param_initialized = true;
            }

            double kappa = (double)Kappa;

            double r;

            if (kappa < 1d) {
                while (true) {
                    double u1 = random.NextUniformOpenInterval01();
                    double u2 = random.NextUniformOpenInterval01();

                    double v = randam_gen_param.c2 * u1;

                    if (v <= 1) {
                        r = randam_gen_param.c1 * double.Pow(v, randam_gen_param.c3);

                        if ((u2 <= (2d - r) / (2d + r)) ||
                            (u2 <= double.Exp(-r))) {

                            break;
                        }
                    }
                    else {
                        r = -double.Log(randam_gen_param.c1 * randam_gen_param.c3 * (randam_gen_param.c2 - v));
                        double y = r / randam_gen_param.c1;

                        if ((u2 * (kappa + y - kappa * y) <= 1d) ||
                            (u2 <= double.Pow(y, kappa - 1d))) {

                            break;
                        }
                    }
                }
            }
            else {
                while (true) {
                    double z = random.NextGaussian();

                    if (randam_gen_param.c2 * z < -1) {
                        continue;
                    }

                    double u = random.NextUniformOpenInterval01();

                    double z2 = z * z;
                    double v = 1d + randam_gen_param.c2 * z;
                    double v3 = v * v * v;

                    if ((u < 1d - 0.0331d * (z2 * z2)) ||
                        (double.Log(u) < (0.5d * z2) + randam_gen_param.c1 * (1d - v3 + double.Log(v3)))) {

                        r = randam_gen_param.c1 * v3;
                        break;
                    }
                }
            }

            double w = r * (double)Theta;

            return w;
        }

        public void SetupRandamGenParams() {
            if ((double)Kappa < 1d) {
                ddouble t, dt, exp_t;

                t = 0.07d + 0.75d * Sqrt(1d - Kappa);

                for (int i = 0; i < 64; i++) {
                    exp_t = Exp(t);
                    dt = (1d - t * (exp_t - 1d) - Kappa) / (1d - (t + 1d) * exp_t);
                    t -= dt;

                    if (Abs(dt) < Abs(t) * 1e-16d) {
                        break;
                    }
                }

                randam_gen_param.c1 = (double)t;
                randam_gen_param.c2 = (double)(1d + Kappa * Exp(-t) / t);
                randam_gen_param.c3 = (double)(1d / Kappa);
            }
            else {
                randam_gen_param.c1 = (double)(Kappa - 1d / 3d);
                randam_gen_param.c2 = (double)(1d / (3d * Sqrt(Kappa - 1d / 3d)));
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Kappa * Theta;

        public override ddouble Median =>
            InverseLowerIncompleteGamma(Kappa, 0.5d) * Theta;

        public override ddouble Mode => Kappa >= 1d ? (Kappa - 1d) * Theta : 0d;

        public override ddouble Variance => Kappa * Theta * Theta;

        public override ddouble Skewness => 2d / Sqrt(Kappa);

        public override ddouble Kurtosis => 6d / Kappa;

        public override ddouble Entropy =>
            Kappa + Log(Theta) + LogGamma(Kappa) - (Kappa - 1d) * Digamma(Kappa);

        public static GammaDistribution operator *(GammaDistribution dist, ddouble k) {
            return new(dist.Kappa, dist.Theta * k);
        }

        public static GammaDistribution operator /(GammaDistribution dist, ddouble k) {
            return new(dist.Kappa, dist.Theta / k);
        }

        public static (GammaDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (GammaDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble kappa = t / (1d - t);

                    try {
                        GammaDistribution dist = new(kappa, 1d);
                        return QuantileScaleFitter<GammaDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble kappa = t / (1d - t);
                GammaDistribution dist = new(kappa, 1d);

                return QuantileScaleFitter<GammaDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(GammaDistribution).Name}[kappa={Kappa},theta={Theta}]";
        }

        public override string Formula => "p(x; kappa, theta) := u^(kappa - 1) * exp(-u) / (gamma(kappa) * theta), u = x / theta";
    }
}
