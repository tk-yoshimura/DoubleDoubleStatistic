using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class NoncentralStudentTDistribution : ContinuousDistribution {

        public ddouble Nu { get; }
        public ddouble Mu { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private const int series_maxiter = 1024;
        private readonly ddouble gc, pdf_norm, pdf_b_scale, nu_inv, nu_half, power;
        private readonly double zero_thr;

        private const int cache_samples = 512;
        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null;

        private readonly ChiSquareDistribution randam_gen_chisq_dist;

        public NoncentralStudentTDistribution(ddouble nu, ddouble mu) {
            ParamAssert.ValidateShape(nameof(nu), ParamAssert.IsFinitePositive(nu));
            ParamAssert.ValidateNonCentricity(nameof(mu), IsFinite(mu));

            Nu = nu;
            Mu = mu;

            ddouble c = Sqrt(nu * PI);

            gc = nu < 70d
                ? Gamma((nu - 1d) / 2d) / Gamma(nu / 2d)
                : Exp(LogGamma((nu - 1d) / 2d) - LogGamma(nu / 2d));

            pdf_b_scale = nu < 70d
                ? Gamma(Nu / 2d + 1d) / Gamma((Nu + 1d) / 2d)
                : Exp(LogGamma(Nu / 2d + 1d) - LogGamma((Nu + 1d) / 2d));

            pdf_norm = nu < 70d
                ? Gamma((nu + 1d) / 2d) / (Gamma(nu / 2d) * c)
                : Exp(LogGamma((nu + 1d) / 2d) - LogGamma(nu / 2d)) / c;

            nu_inv = 1d / nu;
            nu_half = nu / 2d;
            power = -(nu + 1d) / 2d;

            const int zero_thr_log = 710;
            zero_thr = nu < 0.5
                ? double.ScaleB(1, 1000)
                : double.Exp((double)(((nu + 1d) * Log(nu) + 2d * zero_thr_log) / (2d * nu + 2d)));

            randam_gen_chisq_dist = new(nu);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (Abs(x) >= zero_thr) {
                return 0d;
            }

            ddouble u = 1d + x * x * nu_inv;
            ddouble v = Pow(u, power);
            ddouble pdf = pdf_norm * v * PDFScale(x);

            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (Abs(x) >= zero_thr) {
                    return IsNegative(x) ? 0d : 1d;
                }

                ddouble cdf = IsPositive(x) ? CDFKernel(x, Mu) : 1d - CDFKernel(x, -Mu);
                cdf = Clamp(cdf, 0d, 1d);

                return cdf;
            }
            else {
                if (Abs(x) >= zero_thr) {
                    return IsNegative(x) ? 1d : 0d;
                }

                ddouble cdf = IsPositive(x) ? 1d - CDFKernel(x, Mu) : CDFKernel(x, -Mu);
                cdf = Clamp(cdf, 0d, 1d);

                return cdf;
            }
        }

        private ddouble PDFScale(ddouble x) {
            ddouble mux = Mu * x, x2nu = x * x + Nu;
            ddouble u = Square(mux) / (2d * x2nu);

            ddouble da = u * (Nu + 1d);
            ddouble db = u * ((Nu + 2d) / 3d);

            ddouble a = 1d, b = 1d;

            for (int i = 1; i <= series_maxiter; i++) {
                a += da;
                b += db;

                if ((da <= a * 1e-30 || !IsFinite(da)) &&
                    (db <= b * 1e-30 || !IsFinite(db))) {

                    break;
                }
                if (i >= series_maxiter) {
                    throw new ArithmeticException($"{this}: pdf calculation not convergence.");
                }

                da *= u * ((Nu + 1d) * 0.5d + i) / ((i + 1) * (i + 0.5d));
                db *= u * ((Nu + 2d) * 0.5d + i) / ((i + 1) * (i + 1.5d));
            }

            ddouble s = a + mux * Sqrt(2d / x2nu) * pdf_b_scale * b;

            if (s < a * 4e-28) {
                return 0d;
            }

            ddouble y = Exp(Mu * Mu * -0.5d) * s;

            y = Max(0d, y);

            return y;
        }

        private ddouble CDFKernel(ddouble x, ddouble mu) {
            ddouble x2 = x * x, v = x2 / (x2 + Nu);

            ddouble u = mu * mu * 0.5;

            ddouble r1 = 1d, r2 = mu * Sqrt2 / Sqrt(PI);

            ddouble beta1_0 = IncompleteBetaRegularized(v, 0.5d, nu_half);
            ddouble beta1_1 = IncompleteBetaRegularized(v, 1.5d, nu_half);

            ddouble beta2_0 = IncompleteBetaRegularized(v, 1d, nu_half);
            ddouble beta2_1 = IncompleteBetaRegularized(v, 2d, nu_half);

            ddouble s = r1 * beta1_0 + r2 * beta2_0;
            ddouble a1 = 0.5d, a2 = 1d, c1 = nu_half - 0.5d, c2 = nu_half;

            for (int i = 0; i <= series_maxiter; i++) {
                r1 *= u / (i + 1d);
                r2 *= u / (i + 1.5d);

                ddouble ds1 = r1 * beta1_1;
                ddouble ds2 = r2 * beta2_1;

                ddouble ds = ds1 + ds2;

                s += ds;

                if ((Abs(ds1) <= Abs(s) * 1e-30 || !IsFinite(ds1)) &&
                    (Abs(ds2) <= Abs(s) * 1e-30 || !IsFinite(ds2))) {

                    break;
                }
                if (i >= series_maxiter) {
                    throw new ArithmeticException($"{this}: cdf calculation not convergence.");
                }

                a1 += 1d;
                a2 += 1d;
                c1 += 1d;
                c2 += 1d;

                if (beta1_1 > 1e-12) {
                    (beta1_1, beta1_0) = (((a1 + c1 * v) * beta1_1 - c1 * v * beta1_0) / a1, beta1_1);
                }
                else {
                    Debug.WriteLine(
                        "reset recurr incomp.beta: \n" +
                        $"{((a1 + c1 * v) * beta1_1 - c1 * v * beta1_0) / a1} -> {IncompleteBetaRegularized(v, i + 2.5d, nu_half)}"
                    );

                    (beta1_1, beta1_0) = (IncompleteBetaRegularized(v, i + 2.5d, nu_half), beta1_1);
                }

                if (beta2_1 > 1e-12) {
                    (beta2_1, beta2_0) = (((a2 + c2 * v) * beta2_1 - c2 * v * beta2_0) / a2, beta2_1);
                }
                else {
                    Debug.WriteLine(
                        "reset recurr incomp.beta: \n" +
                        $"{((a2 + c2 * v) * beta2_1 - c2 * v * beta2_0) / a2} -> {IncompleteBetaRegularized(v, i + 3d, nu_half)}"
                    );

                    (beta2_1, beta2_0) = (IncompleteBetaRegularized(v, i + 3d, nu_half), beta2_1);
                }
            }

            ddouble y = (Erfc(mu * sqrt2_inv) + Exp(-u) * s) * 0.5d;

            if (1d - y < 1e-29) {
                y = 1d;
            }

            return y;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return interval == Interval.Lower ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return NegativeInfinity;
                }

                this.quantile_lower_builder ??= new QuantileBuilder(-1d, 1d, t => CDF(t / (1d - Abs(t)), Interval.Lower), samples: cache_samples);

                (ddouble t, ddouble t0, ddouble t1) = quantile_lower_builder.Estimate(p);

                ddouble x = t / (1d - Abs(t));
                ddouble x0 = t0 / (1d - Abs(t0));
                ddouble x1 = t1 / (1d - Abs(t1));

                for (int i = 0; i < 64; i++) {
                    ddouble y = CDF(x, Interval.Lower), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    if (IsNegativeInfinity(x0) && IsNegative(x) && IsPositive(dx)) {
                        x = Clamp(x - dx, x * 2d, x1);
                    }
                    else if (IsPositiveInfinity(x1) && IsPositive(x) && IsNegative(dx)) {
                        x = Clamp(x - dx, x0, x * 2d);
                    }
                    else {
                        x = Clamp(x - dx, x0, x1);
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

                this.quantile_upper_builder ??= new QuantileBuilder(1d, -1d, t => CDF(t / (1d - Abs(t)), Interval.Upper), samples: cache_samples);

                (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p);

                ddouble x = t / (1d - Abs(t));
                ddouble x0 = t0 / (1d - Abs(t0));
                ddouble x1 = t1 / (1d - Abs(t1));

                for (int i = 0; i < 64; i++) {
                    ddouble y = CDF(x, Interval.Upper), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    if (IsPositiveInfinity(x0) && IsPositive(x) && IsPositive(dx)) {
                        x = Clamp(x + dx, x1, x * 2d);
                    }
                    else if (IsNegativeInfinity(x1) && IsNegative(x) && IsNegative(dx)) {
                        x = Clamp(x + dx, x * 2d, x0);
                    }
                    else {
                        x = Clamp(x + dx, x1, x0);
                    }

                    if (Abs(dx) <= Abs(x) * 1e-29) {
                        break;
                    }
                }

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextGaussian();
            double v = randam_gen_chisq_dist.Sample(random);

            double r = (u + (double)Mu) / double.Sqrt(v * (double)nu_inv);

            return r;
        }

        public override bool Symmetric => Mu == 0d;

        public override ddouble Mean => Nu > 1d
            ? Mu * Sqrt(Nu * 0.5d) * gc
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode =>
            mode ??= BisectionMinimizeSearch1D.Search(t => -PDF(t), (Mu * Sqrt(Nu / (Nu + 2.5d)), Mu * Sqrt(Nu / (Nu + 1d))), iter: 1024);

        public override ddouble Variance => Nu > 2d
            ? Nu * (1d + Mu * Mu) / (Nu - 2d) - Mu * Mu * Nu * 0.5d * Square(gc)
            : NaN;

        public override ddouble Skewness {
            get {
                if (Nu <= 3d) {
                    return NaN;
                }

                ddouble variance = Variance;

                return Mu * Sqrt(Nu * 0.5d) * gc * (Nu * (Mu * Mu + 2d * Nu - 3d) / ((Nu - 3d) * (Nu - 2d)) - 2d * variance) / ExMath.Pow3d2(variance);
            }
        }

        public override ddouble Kurtosis {
            get {
                if (Nu <= 4d) {
                    return NaN;
                }

                ddouble variance = Variance;
                ddouble mu2 = Mu * Mu;

                return (Nu * Nu * (mu2 * mu2 + 6d * mu2 + 3d) / ((Nu - 4d) * (Nu - 2d)) -
                    mu2 * Nu * Square(gc) * 0.5d * (Nu * (mu2 * (Nu + 1d) + (9d * Nu - 15d)) / ((Nu - 3d) * (Nu - 2d)) - 3d * variance)) / Square(variance) - 3d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(NoncentralStudentTDistribution).Name}[nu={Nu},mu={Mu}]";
        }

        public override string Formula =>
            "p(x; nu, mu) := (1 + x^2 / nu)^(- (nu + 1) / 2) * gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(pi * nu)) * exp(-mu^2 / 2) * (A_nu(x; mu) + B_nu(x; mu)), " +
            "A_nu(x; mu) := Hypergeometric1F1((nu + 1) / 2; 1 / 2; mu^2 * x^2 / (2 * (x^2 + nu))), " +
            "B_nu(x; mu) := (sqrt(2) * mu * x / sqrt(x^2 + nu)) * gamma(nu / 2 + 1) / gamma((nu + 1) / 2) * Hypergeometric1F1(nu / 2 + 1; 3 / 2; mu^2 * x^2 / (2 * (x^2 + nu)))";
    }
}
