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
    public class HotellingTSquareDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<HotellingTSquareDistribution> {

        public ddouble P { get; }
        public ddouble M { get; }

        private readonly ddouble scale, scale_inv;
        private readonly SnedecorFDistribution fdist;

        public HotellingTSquareDistribution(ddouble p, ddouble m) {
            ParamAssert.ValidateShape(nameof(p), ParamAssert.IsFinitePositive(p));
            ParamAssert.ValidateShape(nameof(m), ParamAssert.IsFinitePositive(m));

            P = p;
            M = m;

            fdist = new SnedecorFDistribution(p, 1d - p + m);
            scale = (m * p) / (m - p + 1d);
            scale_inv = (m - p + 1d) / (m * p);
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = fdist.PDF(x * scale_inv) * scale_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble cdf = fdist.CDF(x * scale_inv, interval);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            ddouble x = fdist.Quantile(p, interval) * scale;

            return x;
        }

        public override double Sample(Random random) {
            double v = fdist.Sample(random) * (double)scale;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => fdist.Mean * scale;

        public override ddouble Median => fdist.Median * scale;

        public override ddouble Mode => fdist.Mode * scale;

        public override ddouble Variance => fdist.Variance * scale * scale;

        public override ddouble Skewness => fdist.Skewness;

        public override ddouble Kurtosis => fdist.Kurtosis;

        public override ddouble Entropy => fdist.Entropy + Log(scale);

        public static (HotellingTSquareDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (HotellingTSquareDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = BisectionMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble p = t.u / (1d - t.u);
                    ddouble m = t.v / (1d - t.v);

                    try {
                        HotellingTSquareDistribution dist = new(p, m);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-4d, 1e-4d), (1000d / 1001d, 1000d / 1001d)), iter: 64
            );

            try {
                ddouble p = u / (1d - u);
                ddouble m = v / (1d - v);
                HotellingTSquareDistribution dist = new(p, m);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(HotellingTSquareDistribution).Name}[p={P},m={M}]";
        }

        public override string Formula => "p(x; p, m) := snedecor_f(x * (m - p + 1) / (m * p); p, 1 - p + m)";
    }
}
