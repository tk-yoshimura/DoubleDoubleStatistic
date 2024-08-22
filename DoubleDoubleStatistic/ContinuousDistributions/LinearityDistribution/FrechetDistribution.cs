﻿using DoubleDouble;
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
    public class FrechetDistribution : LinearityDistribution<FrechetDistribution>,
        IAdditionOperators<FrechetDistribution, ddouble, FrechetDistribution>,
        ISubtractionOperators<FrechetDistribution, ddouble, FrechetDistribution>,
        IMultiplyOperators<FrechetDistribution, ddouble, FrechetDistribution>,
        IDivisionOperators<FrechetDistribution, ddouble, FrechetDistribution>,
        IFittableContinuousDistribution<FrechetDistribution> {

        public ddouble Alpha { get; }
        public ddouble Mu { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_norm, alpha_inv, theta_inv;

        public FrechetDistribution(ddouble alpha) : this(alpha: alpha, mu: 0d, theta: 1d) { }

        public FrechetDistribution(ddouble alpha, ddouble theta) : this(alpha: alpha, mu: 0d, theta: theta) { }

        public FrechetDistribution(ddouble alpha, ddouble mu, ddouble theta) {
            ParamAssert.ValidateShape(nameof(alpha), ParamAssert.IsFinitePositive(alpha));
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(theta), ParamAssert.IsFinitePositive(theta));

            Alpha = alpha;
            Mu = mu;
            Theta = theta;

            pdf_norm = alpha / theta;
            alpha_inv = 1d / alpha;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * theta_inv;

            if (u <= 0d) {
                return 0d;
            }
            if (IsNaN(u)) {
                return NaN;
            }

            ddouble v = Log2(u) * Alpha;

            ddouble pdf = pdf_norm * Pow2(-Pow2(-v) * LbE - v) / u;

            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }

                ddouble v = Pow(u, Alpha);

                if (!IsFinite(v)) {
                    return 1d;
                }

                ddouble cdf = Exp(-1d / v);

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }

                ddouble v = Pow(u, Alpha);

                if (!IsFinite(v)) {
                    return 0d;
                }

                ddouble cdf = -Expm1(-1d / v);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p <= 0d) {
                    return Mu;
                }
                if (p >= 1d) {
                    return PositiveInfinity;
                }

                ddouble v = 1d / -Log(p);

                if (IsNegative(v)) {
                    return PositiveInfinity;
                }

                ddouble u = Pow(v, alpha_inv);
                ddouble x = Mu + u * Theta;

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }
                if (p >= 1d) {
                    return Mu;
                }

                ddouble v = 1d / -Log1p(-p);

                if (IsNegative(v)) {
                    return 0d;
                }

                ddouble u = Pow(v, alpha_inv);
                ddouble x = Mu + u * Theta;

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(1d / -double.Log(u), (double)alpha_inv);
            double w = (double)Mu + v * (double)Theta;

            return w;
        }

        public override (ddouble min, ddouble max) Support => (Mu, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d)
            ? Mu + Theta * Gamma(1d - alpha_inv)
            : PositiveInfinity;

        public override ddouble Median =>
            Mu + Theta / Pow(Ln2, alpha_inv);

        public override ddouble Mode =>
            Mu + Pow1p(-1d / (1d + Alpha), alpha_inv) * Theta;

        public override ddouble Variance => (Alpha > 2d)
            ? Square(Theta) * (Gamma(1d - 2d * alpha_inv) - Square(Gamma(1d - alpha_inv)))
            : PositiveInfinity;

        public override ddouble Skewness {
            get {
                if (Alpha <= 3d) {
                    return PositiveInfinity;
                }

                ddouble g1 = Gamma(1d - alpha_inv), g2 = Gamma(1d - 2d * alpha_inv), g3 = Gamma(1d - 3d * alpha_inv);

                return (g3 - 3d * g2 * g1 + 2d * Cube(g1)) / ExMath.Pow3d2(g2 - Square(g1));
            }
        }

        public override ddouble Kurtosis {
            get {
                if (Alpha <= 4d) {
                    return PositiveInfinity;
                }

                ddouble g1 = Gamma(1d - alpha_inv), g2 = Gamma(1d - 2d * alpha_inv);
                ddouble g3 = Gamma(1d - 3d * alpha_inv), g4 = Gamma(1d - 4d * alpha_inv);

                return (g4 - 4d * g3 * g1 + 3d * Square(g2)) / Square(g2 - Square(g1)) - 6d;
            }
        }

        public override ddouble Entropy =>
            1d + EulerGamma * (Alpha + 1d) * alpha_inv + Log(Theta * alpha_inv);

        public static FrechetDistribution operator +(FrechetDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Theta);
        }

        public static FrechetDistribution operator -(FrechetDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Theta);
        }

        public static FrechetDistribution operator *(FrechetDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Theta * k);
        }

        public static FrechetDistribution operator /(FrechetDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu / k, dist.Theta / k);
        }

        public static (FrechetDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (FrechetDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble alpha = t / (1d - t);

                    try {
                        FrechetDistribution dist = new(alpha);
                        return QuantileLinearFitter<FrechetDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble alpha = t / (1d - t);
                FrechetDistribution dist = new(alpha);

                return QuantileLinearFitter<FrechetDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(FrechetDistribution).Name}[alpha={Alpha},mu={Mu},theta={Theta}]";
        }

        public override string Formula => "p(x; alpha, mu, sigma) := exp(- u ^ (-alpha)) * alpha / u^(1 + alpha) / sigma, u = (x - mu) / sigma";
    }
}
