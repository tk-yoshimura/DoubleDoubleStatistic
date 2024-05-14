using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;

namespace DoubleDoubleStatistic.Utils {
    public static class QuantileLinearFitter<Distribution>
        where Distribution : LinearityDistribution<Distribution>,
            IMultiplyOperators<Distribution, ddouble, Distribution>,
            IDivisionOperators<Distribution, ddouble, Distribution>,
            IAdditionOperators<Distribution, ddouble, Distribution>,
            ISubtractionOperators<Distribution, ddouble, Distribution> {

        public static (Distribution? dist, ddouble error) Fit(
            Distribution dist_base,
            IEnumerable<double> samples,
            (double min, double max) fitting_quantile_range,
            int quantile_partitions) {

            return Fit(dist_base, samples.Select(v => (ddouble)v).ToArray(), fitting_quantile_range, quantile_partitions);
        }

        public static (Distribution? dist, ddouble error) Fit(
            Distribution dist_base,
            ReadOnlyCollection<double> samples,
            (double min, double max) fitting_quantile_range,
            int quantile_partitions) {

            return Fit(dist_base, samples.Select(v => (ddouble)v).ToArray(), fitting_quantile_range, quantile_partitions);
        }

        public static (Distribution? dist, ddouble error) Fit(
            Distribution dist_base,
            IEnumerable<ddouble> samples,
            (ddouble min, ddouble max) fitting_quantile_range,
            int quantile_partitions) {

            return FitForSortedSamples(dist_base, new ReadOnlyCollection<ddouble>(samples.Sort().ToArray()), fitting_quantile_range, quantile_partitions);
        }

        public static (Distribution? dist, ddouble error) Fit(
            Distribution dist_base,
            ReadOnlyCollection<ddouble> samples,
            (ddouble min, ddouble max) fitting_quantile_range,
            int quantile_partitions) {

            return FitForSortedSamples(dist_base, new ReadOnlyCollection<ddouble>(samples.Sort().ToArray()), fitting_quantile_range, quantile_partitions);
        }

        internal static (Distribution? dist, ddouble error) FitForSortedSamples(
            Distribution dist_base,
            ReadOnlyCollection<ddouble> sorted_samples,
            (ddouble min, ddouble max) fitting_quantile_range,
            int quantile_partitions) {

            if (!(fitting_quantile_range.min < fitting_quantile_range.max && fitting_quantile_range.min >= 0d && fitting_quantile_range.max <= 1d)) {
                throw new ArgumentException("Invalid range: min < max", nameof(fitting_quantile_range));
            }

            ArgumentOutOfRangeException.ThrowIfLessThan(quantile_partitions, 10, nameof(quantile_partitions));

            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = sorted_samples.SortedQuantile(qs).ToArray();

            return FitForQuantiles(dist_base, qs, ys);
        }

        internal static (Distribution? dist, ddouble error) FitForQuantiles(Distribution dist_base, ddouble[] qs, ddouble[] ys) {

            Debug.Assert(qs.Length == ys.Length);

            List<ddouble> xs = qs.Select(q => dist_base.Quantile(q)).ToList();

            int n = 0;
            ddouble sum_x = 0d, sum_y = 0d, sum_x2 = 0d, sum_xy = 0d;

            for (int i = 0; i < qs.Length; i++) {
                ddouble x = xs[i], y = ys[i];

                n++;
                sum_x += x;
                sum_y += y;
                sum_x2 += x * x;
                sum_xy += x * y;
            }

            try {
                ddouble r = 1d / (sum_x * sum_x - n * sum_x2);
                ddouble scale = (sum_x * sum_y - n * sum_xy) * r;
                ddouble mu = (sum_x * sum_xy - sum_x2 * sum_y) * r;

                Distribution dist = dist_base * scale + mu;

                ddouble error = xs.Select((x, idx) => ddouble.Square(x * scale + mu - ys[idx])).Mean();

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, ddouble.NaN);
            }
        }
    }
}
