﻿using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class NoncentralSnedecorFDistribution : ContinuousDistribution {

        public ddouble N { get; }
        public ddouble M { get; }
        public ddouble Lambda { get; }

        private const int series_maxiter = 8192;
        private const int cache_samples = 512;
        private QuantileBuilder? quantile_lower_builder = null, quantile_upper_builder = null;

        private readonly ChiSquareDistribution randam_gen_chisq_dist;
        private readonly NoncentralChiSquareDistribution randam_gen_noncchisq_dist;

        private ddouble cdf_threshold = 0d;

        public NoncentralSnedecorFDistribution(ddouble n, ddouble m, ddouble lambda) {
            ParamAssert.ValidateShape(nameof(n), ParamAssert.IsFinitePositive(n));
            ParamAssert.ValidateShape(nameof(m), ParamAssert.IsFinitePositive(m));
            ParamAssert.ValidateNonCentricity(nameof(lambda), ParamAssert.IsFinitePositive(lambda));

            N = n;
            M = m;
            Lambda = lambda;

            randam_gen_chisq_dist = new(m);
            randam_gen_noncchisq_dist = new(n, lambda);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble a = N / M, b = M / (N * x + M), r = a * b * x * Lambda * 0.5d;
            ddouble v = Pow(a, N * 0.5d) * Pow(b, (N + M) * 0.5d) * Pow(x, N * 0.5d - 1);
            ddouble beta = Beta(M * 0.5d, N * 0.5);

            ddouble exp = Lambda * LbE * -0.5d;

            ddouble u = 0d;

            for (int i = 0; i <= series_maxiter; i++) {
                ddouble du = v / beta;

                u += du;

                if (Abs(du) <= Abs(u) * 1e-30 || !IsFinite(du)) {
                    break;
                }

                if (i >= series_maxiter) {
                    throw new ArithmeticException($"{this}: pdf calculation not convergence.");
                }

                v *= r / (i + 1);
                beta *= (N * 0.5 + i) / ((M + N) * 0.5d + i);

                // Overflow avoidance by rescaling
                if (double.ILogB((double)u) >= 16 || double.ILogB((double)v) >= 16) {
                    (u, v) = (Ldexp(u, -16), Ldexp(v, -16));
                    exp += 16d;
                }
            }

            ddouble pdf = Pow2(exp) * u;
            pdf = IsFinite(pdf) ? Max(pdf, 0d) : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d) {
                return interval == Interval.Lower ? 0d : 1d;
            }

            if (IsPositiveInfinity(x)) {
                return interval == Interval.Lower ? 1d : 0d;
            }

            if (cdf_threshold == 0d) {
                for (cdf_threshold = Ldexp(1, -64); double.ILogB((double)cdf_threshold) < 64; cdf_threshold *= 2d) {
                    if (CDFLower(cdf_threshold) > 0.25) {
                        break;
                    }
                }
            }

            if (x < cdf_threshold) {
                ddouble cdf = CDFLower(x);

                if (interval == Interval.Upper) {
                    cdf = 1d - cdf;
                }

                return cdf;
            }
            else {
                ddouble cdf = CDFUpper(x);

                if (interval == Interval.Lower) {
                    cdf = 1d - cdf;
                }

                return cdf;
            }
        }


        private ddouble CDFLower(ddouble x) {
            ddouble lambda_half = Lambda * 0.5d, n_half = N * 0.5d, m_half = M * 0.5d;
            ddouble exp = -lambda_half * LbE;
            ddouble r = 1d;
            ddouble u = N * x / (N * x + M);

            ddouble beta0 = IncompleteBetaRegularized(u, n_half, m_half);
            ddouble beta1 = IncompleteBetaRegularized(u, n_half + 1d, m_half);

            ddouble s = beta0;
            ddouble a = n_half, c = a + m_half - 1d;

            for (int i = 1; i <= series_maxiter; i++) {
                r *= lambda_half / i;
                ddouble ds = r * beta1;

                s += ds;

                if (Abs(ds) <= Abs(s) * 1e-30 || !IsFinite(ds)) {
                    break;
                }

                if (i >= series_maxiter) {
                    throw new ArithmeticException($"{this}: cdf calculation not convergence.");
                }

                a += 1d;
                c += 1d;

                if (beta1 > 1e-12) {
                    (beta1, beta0) = (((a + c * u) * beta1 - c * u * beta0) / a, beta1);
                }
                else {
                    Debug.WriteLine(
                        "reset recurr incomp.beta: \n" +
                        $"{((a + c * u) * beta1 - c * u * beta0) / a} -> {IncompleteBetaRegularized(u, n_half + (i + 1), m_half)}"
                    );

                    (beta1, beta0) = (IncompleteBetaRegularized(u, n_half + (i + 1), m_half), beta1);
                }

                // Overflow avoidance by rescaling
                if (double.ILogB((double)s) >= 16 || double.ILogB((double)r) >= 16) {
                    (s, r) = (Ldexp(s, -16), Ldexp(r, -16));
                    exp += 16d;
                }
            }

            ddouble cdf = Min(1d, Pow2(exp) * s);
            return cdf;
        }
        private ddouble CDFUpper(ddouble x) {
            ddouble lambda_half = Lambda * 0.5d, n_half = N * 0.5d, m_half = M * 0.5d;
            ddouble exp = -lambda_half * LbE;
            ddouble r = 1d;
            ddouble u = M / (N * x + M), v = N * x / (N * x + M);

            ddouble beta0 = IncompleteBetaRegularized(u, m_half, n_half);
            ddouble beta1 = IncompleteBetaRegularized(u, m_half, n_half + 1d);

            ddouble s = beta0;
            ddouble b = n_half, c = m_half + b - 1d;

            for (int i = 1; i <= series_maxiter; i++) {
                r *= lambda_half / i;
                ddouble ds = r * beta1;

                s += ds;

                if (Abs(ds) <= Abs(s) * 1e-30 || !IsFinite(ds)) {
                    break;
                }

                if (i >= series_maxiter) {
                    throw new ArithmeticException($"{this}: cdf calculation not convergence.");
                }

                b += 1d;
                c += 1d;

                if (beta1 > 1e-12) {
                    (beta1, beta0) = (((b + c * v) * beta1 - c * v * beta0) / b, beta1);
                }
                else {
                    Debug.WriteLine(
                        "reset recurr incomp.beta: \n" +
                        $"{((b + c * v) * beta1 - c * v * beta0) / b} -> {IncompleteBetaRegularized(u, m_half, n_half + (i + 1))}"
                    );

                    (beta1, beta0) = (IncompleteBetaRegularized(u, m_half, n_half + (i + 1)), beta1);
                }

                // Overflow avoidance by rescaling
                if (double.ILogB((double)s) >= 16 || double.ILogB((double)r) >= 16) {
                    (s, r) = (Ldexp(s, -16), Ldexp(r, -16));
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
                    1d, 0d,
                    t => CDF((1d - t) / t, Interval.Lower),
                    cache_samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_lower_builder.Estimate(p);

                ddouble x = (1d - t) / t;
                ddouble x0 = (1d - t0) / t0;
                ddouble x1 = (1d - t1) / t1;

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
                    0d, 1d,
                    t => CDF((1d - t) / t, Interval.Upper),
                    cache_samples
                );

                (ddouble t, ddouble t0, ddouble t1) = quantile_upper_builder.Estimate(p);

                ddouble x = (1d - t) / t;
                ddouble x0 = (1d - t0) / t0;
                ddouble x1 = (1d - t1) / t1;

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
            double t = randam_gen_noncchisq_dist.Sample(random) * (double)M;
            double s = randam_gen_chisq_dist.Sample(random) * (double)N;

            double r = t / s;

            return r;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => (M > 2d)
            ? (M * (Lambda + N)) / ((M - 2d) * N)
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble t = BisectionMinimizeSearch1D.Search(t => -PDF(t / (1d - t)), (0.125d, 1d), iter: 1024);
                ddouble xp = t / (1d - t);

                ddouble xn = BisectionMinimizeSearch1D.Search(t => -PDF(t), (0d, 0.15d), iter: 1024);

                ddouble vp = PDF(xp), vn = PDF(xn);

                mode ??= (vp >= vn) ? xp : xn;

                return mode.Value;
            }
        }

        public override ddouble Variance => (M > 4d)
            ? 2d * M * M * ((M - 2d) * (2d * Lambda + N) + Square(Lambda + N)) / ((M - 4d) * Square(N * (M - 2d)))
            : NaN;

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??= (M > 6d)
            ? IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??= (M > 8d)
            ? IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048)
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 32768);

        public override string ToString() {
            return $"{typeof(NoncentralSnedecorFDistribution).Name}[n={N},m={M},lambda={Lambda}]";
        }

        public override string Formula =>
            "p(x; n, m, lambda) := exp(-lambda / 2) * sum((lambda / 2)^k / (beta(m / 2, n / 2 + k) * k!) * (n / m)^(n / 2 + k) * (m / (n * x + m))^((n + m) / 2 + k) * x^(n / 2 - 1 + k), k, 0, inf)";
    }
}
