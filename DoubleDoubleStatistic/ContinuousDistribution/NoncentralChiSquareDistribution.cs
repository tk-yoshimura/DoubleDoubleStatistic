using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Utils;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class NoncentralChiSquareDistribution : ContinuousDistribution {

        public ddouble Nu { get; }
        public ddouble Lambda { get; }

        private readonly ddouble c;

        private const int series_maxiter = 1024;
        private const int cache_samples = 512;
        private QuantileBuilder quantile_lower_builder = null, quantile_upper_builder = null;

        public NoncentralChiSquareDistribution(ddouble nu, ddouble lambda) {
            ValidateShape(nu, nu => nu > 0d);
            ValidateShape(lambda, lambda => lambda > 0d);

            Nu = nu;
            Lambda = lambda;

            c = nu * 0.5d - 1d;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble u = Max(Epsilon, x / Lambda);

            ddouble pdf = Exp((-Square(Sqrt(x) - Sqrt(Lambda)) + Log(u) * c) * 0.5d) * BesselI(c, Sqrt(x * Lambda), scale: true) * 0.5d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d) {
                return interval == Interval.Lower ? 0d : 1d;
            }

            if (IsPositiveInfinity(x) || (x > Mean && PDF(x) == 0d)) {
                return interval == Interval.Lower ? 1d : 0d;
            }

            ddouble nu_half = Nu * 0.5d, lambda_half = Lambda * 0.5d, x_half = x * 0.5d;
            ddouble f = 1d, u = nu_half, v = Gamma(u);
            ddouble lnx_half = Log(x_half);

            if (x < Mean) {
                ddouble g = LowerIncompleteGammaRegularized(u, x_half);
                ddouble s = g;

                for (int i = 1; i <= series_maxiter; i++) {
                    f *= lambda_half / i;
                    g -= Exp(lnx_half * u - x_half) / (u * v);

                    ddouble ds = f * g;

                    s += ds;

                    if (ds <= s * 1e-30 || !IsFinite(ds)) {
                        break;
                    }

                    if (i >= series_maxiter) {
                        throw new ArithmeticException($"{this}: cdf calculation not convergence.");
                    }

                    v *= u;
                    u += 1d;
                }

                ddouble cdf = Min(1d, Exp(-lambda_half) * s);

                if (interval == Interval.Upper) {
                    cdf = 1d - cdf;
                }

                return cdf;
            }
            else {
                ddouble g = UpperIncompleteGammaRegularized(u, x_half);
                ddouble s = g;

                for (int i = 1; i <= series_maxiter; i++) {
                    f *= lambda_half / i;
                    g += Exp(lnx_half * u - x_half) / (u * v);

                    ddouble ds = f * g;

                    s += ds;

                    if (ds <= s * 1e-30 || !IsFinite(ds)) {
                        break;
                    }

                    if (i >= series_maxiter) {
                        throw new ArithmeticException($"{this}: cdf calculation not convergence.");
                    }

                    v *= u;
                    u += 1d;
                }

                ddouble cdf = Min(1d, Exp(-lambda_half) * s);

                if (interval == Interval.Lower) {
                    cdf = 1d - cdf;
                }

                return cdf;
            }
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

                    if (Abs(dx / x) < 1e-29 || Abs(dx) < Epsilon) {
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

                    if (Abs(dx / x) < 1e-29 || Abs(dx) < Epsilon) {
                        break;
                    }
                }

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble t = MaximumFinder.BisectionFind(0d, 0.875d, t => PDF((1d - t) / t));
                ddouble xp = (1d - t) / t;
                ddouble xn = MaximumFinder.BisectionFind(0d, 0.15d, PDF);

                ddouble vp = PDF(xp), vn = PDF(xn);

                mode ??= (vp >= vn) ? xp : xn;

                return mode.Value;
            }
        }

        public override ddouble Mean => Nu + Lambda;

        public override ddouble Variance => 2d * Nu + 4d * Lambda;

        public override ddouble Skewness => (Nu + 3d * Lambda) * ExMath.Pow3d2(2d / (Nu + 2 * Lambda));

        public override ddouble Kurtosis => 12d * (Nu + 4d * Lambda) / Square(Nu + 2d * Lambda);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 65536);

        public override string ToString() {
            return $"{typeof(NoncentralChiSquareDistribution).Name}[nu={Nu},lambda={Lambda}]";
        }

        public override string Formula => "p(x; nu, lambda) := (x / lambda)^(nu / 4 - 1 / 2) * exp(-(x + lambda) / 2) * bessel_i(nu / 2 - 1, sqrt(x * lambda)) / 2";
    }
}
