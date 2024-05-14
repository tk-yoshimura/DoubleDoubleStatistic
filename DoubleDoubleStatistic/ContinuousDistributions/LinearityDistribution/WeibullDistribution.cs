using DoubleDouble;
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
    public class WeibullDistribution : LinearityDistribution<WeibullDistribution>,
        IAdditionOperators<WeibullDistribution, ddouble, WeibullDistribution>,
        ISubtractionOperators<WeibullDistribution, ddouble, WeibullDistribution>,
        IMultiplyOperators<WeibullDistribution, ddouble, WeibullDistribution>,
        IDivisionOperators<WeibullDistribution, ddouble, WeibullDistribution>,
        IFittableDistribution<WeibullDistribution> {

        public ddouble Alpha { get; }
        public ddouble Mu { get; }
        public ddouble Theta { get; }

        private readonly ddouble pdf_norm, alpha_inv, theta_inv;

        public WeibullDistribution(ddouble alpha) : this(alpha: alpha, mu: 0d, theta: 1d) { }

        public WeibullDistribution(ddouble alpha, ddouble theta) : this(alpha: alpha, mu: 0d, theta: theta) { }

        public WeibullDistribution(ddouble alpha, ddouble mu, ddouble theta) {
            ValidateShape(alpha, alpha => alpha > 0d);
            ValidateLocation(mu);
            ValidateScale(theta);

            Alpha = alpha;
            Mu = mu;
            Theta = theta;

            pdf_norm = alpha / theta;
            alpha_inv = 1d / alpha;
            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            if (u <= 0d) {
                return Alpha < 1d ? PositiveInfinity : Alpha == 1d ? theta_inv : 0d;
            }
            if (IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble v = Log2(u) * Alpha;

            ddouble pdf = pdf_norm * Pow2(-Pow2(v) * LbE + v) / u;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * theta_inv;

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(u)) {
                    return 1d;
                }

                ddouble cdf = -Expm1(-Pow(u, Alpha));

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = Exp(-Pow(u, Alpha));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p == 1d) {
                    return PositiveInfinity;
                }

                ddouble u = Pow(-Log1p(-p), alpha_inv);
                ddouble x = Mu + u * Theta;

                return x;
            }
            else {
                if (p == 0d) {
                    return PositiveInfinity;
                }

                ddouble u = Pow(-Log(p), alpha_inv);
                ddouble x = Mu + u * Theta;

                if (IsNegative(x)) {
                    return 0d;
                }

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = double.Pow(-double.Log(u), (double)alpha_inv);
            double w = (double)Mu + v * (double)Theta;

            return w;
        }

        public override (ddouble min, ddouble max) Support => (Mu, PositiveInfinity);

        public override ddouble Mean =>
            Mu + Theta * Gamma(1d + alpha_inv);

        public override ddouble Median =>
            Mu + Theta * Pow(Ln2, alpha_inv);

        public override ddouble Mode => Alpha <= 1d
            ? Mu
            : Mu + Theta * Pow((Alpha - 1d) * alpha_inv, alpha_inv);

        public override ddouble Variance =>
            Theta * Theta * (Gamma(1d + 2d * alpha_inv) - Square(Gamma(1d + alpha_inv)));

        public override ddouble Skewness {
            get {
                ddouble mu = Gamma(1d + alpha_inv), var = Gamma(1d + 2d * alpha_inv) - Square(mu);

                return (Gamma(1d + 3d * alpha_inv) - 3d * mu * var - Cube(mu)) / ExMath.Pow3d2(var);
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mu = Gamma(1d + alpha_inv), var = Gamma(1d + 2d * alpha_inv) - Square(mu);

                return (Gamma(1d + 4d * alpha_inv)
                    - 4d * mu * (Gamma(1d + 3d * alpha_inv) - 3d * mu * var - Cube(mu))
                    - 6d * Square(mu) * var
                    - Square(Square(mu))) /
                    Square(var) - 3d;
            }
        }

        public override ddouble Entropy => 1d + EulerGamma * (1d - alpha_inv) + Log(Theta * alpha_inv);

        public static WeibullDistribution operator +(WeibullDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Theta);
        }

        public static WeibullDistribution operator -(WeibullDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Theta);
        }

        public static WeibullDistribution operator *(WeibullDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Theta * k);
        }

        public static WeibullDistribution operator /(WeibullDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu / k, dist.Theta / k);
        }

        public static (WeibullDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (WeibullDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble alpha = t / (1d - t);

                    try {
                        WeibullDistribution dist = new(alpha);
                        return QuantileLinearFitter<WeibullDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            ddouble alpha = t / (1d - t);
            WeibullDistribution dist = new(alpha);

            return QuantileLinearFitter<WeibullDistribution>.FitForQuantiles(dist, qs, ys);
        }

        public override string ToString() {
            return $"{typeof(WeibullDistribution).Name}[alpha={Alpha},mu={Mu},theta={Theta}]";
        }

        public override string Formula => "p(x; alpha, mu, theta) := u^(k - 1) * exp(-u^alpha) * alpha / theta, u = (x - mu) / theta";
    }
}
