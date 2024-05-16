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
    public class JohnsonSBDistribution : LinearityDistribution<JohnsonSBDistribution>,
        IMultiplyOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        IAdditionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        ISubtractionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        IDivisionOperators<JohnsonSBDistribution, ddouble, JohnsonSBDistribution>,
        IFittableContinuousDistribution<JohnsonSBDistribution> {

        public ddouble Gamma { get; }
        public ddouble Delta { get; }
        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private readonly ddouble pdf_norm, sigma_inv;

        public JohnsonSBDistribution(ddouble gamma, ddouble delta) : this(gamma, delta, mu: 0d, sigma: 1d) { }

        public JohnsonSBDistribution(ddouble gamma, ddouble delta, ddouble sigma) : this(gamma, delta, mu: 0d, sigma: sigma) { }

        public JohnsonSBDistribution(ddouble gamma, ddouble delta, ddouble mu, ddouble sigma) {
            ParamAssert.ValidateShape(nameof(gamma), IsFinite(gamma));
            ParamAssert.ValidateShape(nameof(delta), IsFinite(delta));
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Gamma = gamma;
            Delta = delta;
            Mu = mu;
            Sigma = sigma;

            pdf_norm = delta / (sigma * Sqrt(2d * PI));
            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (u <= 0d || u >= 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm / (u * (1d - u)) * Exp(-0.5d * Square(Gamma + Delta * Log(u / (1d - u))));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (u <= 0d) {
                return (interval == Interval.Lower) ? 0d : 1d;
            }
            if (u >= 1d) {
                return (interval == Interval.Lower) ? 1d : 0d;
            }

            ddouble v = Gamma + Delta * Log(u / (1d - u));

            if (interval == Interval.Lower) {
                ddouble cdf = Erfc(-v * sqrt2_inv) * 0.5d;

                return cdf;
            }
            else {
                ddouble cdf = Erfc(v * sqrt2_inv) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble v = interval == Interval.Lower
                ? -InverseErfc(p * 2d) * Sqrt2
                : InverseErfc(p * 2d) * Sqrt2;

            ddouble u = 1d / (Exp((Gamma - v) / Delta) + 1d);
            ddouble x = Mu + Sigma * u;

            return x;
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();
            double v = ErrorFunction.InverseErfc(u * 2d) * double.Sqrt(2d);
            double w = (double)Mu + (double)Sigma / (double.Exp(((double)Gamma - v) / (double)Delta) + 1d);

            return w;
        }

        public override (ddouble min, ddouble max) Support => (Mu, Mu + Sigma);

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble u = 0.5d;

                for (int i = 0; i < 256; i++) {
                    ddouble v = u * (u - 1d);

                    ddouble du = (2d * u - 1d - Delta * (Gamma + Delta * Log(u / (1d - u)))) * v
                        / (2d * v + Delta * Delta);

                    if (!IsFinite(du)) {
                        break;
                    }

                    u = Clamp(u - du, u / 16d, (u + 15d) / 16d);

                    if (Abs(du) <= Abs(u) * 1e-29) {
                        break;
                    }
                }

                ddouble x = mode ??= Mu + u * Sigma;

                return x;
            }
        }

        public override ddouble Median => Quantile(0.5d);

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
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static JohnsonSBDistribution operator *(JohnsonSBDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Delta, dist.Mu * k, dist.Sigma * k);
        }

        public static JohnsonSBDistribution operator /(JohnsonSBDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Delta, dist.Mu / k, dist.Sigma / k);
        }

        public static JohnsonSBDistribution operator +(JohnsonSBDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Delta, dist.Mu + s, dist.Sigma);
        }

        public static JohnsonSBDistribution operator -(JohnsonSBDistribution dist, ddouble s) {
            return new(dist.Gamma, dist.Delta, dist.Mu - s, dist.Sigma);
        }

        public static (JohnsonSBDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (JohnsonSBDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = GridMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble gamma = t.u / (1d - t.u);
                    ddouble delta = t.v / (1d - t.v);

                    try {
                        JohnsonSBDistribution dist = new(gamma, delta);
                        return QuantileLinearFitter<JohnsonSBDistribution>.FitForQuantiles(dist, qs, ys).error;
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-10d, 1e-10d), (10d / 11d, 10d / 11d)), iter: 64
            );

            try {
                ddouble gamma = u / (1d - u);
                ddouble delta = v / (1d - v);
                JohnsonSBDistribution dist = new(gamma, delta);

                return QuantileLinearFitter<JohnsonSBDistribution>.FitForQuantiles(dist, qs, ys);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(JohnsonSBDistribution).Name}[gamma={Gamma},delta={Delta},mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; gamma, delta, mu, sigma) := (delta * exp(-(gamma + delta * log(u / (1 - u)))^2 / 2)) / (sqrt(2 * pi) * u * (1 - u)) / sigma, u = (x - mu) / sigma";
    }
}
