using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using System.Collections.ObjectModel;
using System.Numerics;

namespace DoubleDoubleStatistic.Utils {
    public static class QuantileLinearFitter<Distribution>
        where Distribution : LinearityDistribution<Distribution>,
            IMultiplyOperators<Distribution, ddouble, Distribution>,
            IDivisionOperators<Distribution, ddouble, Distribution>,
            IAdditionOperators<Distribution, ddouble, Distribution>,
            ISubtractionOperators<Distribution, ddouble, Distribution> {

        public static Distribution Fit(
            Distribution dist_base,
            IEnumerable<double> samples,
            (double min, double max) fitting_quantile_range) {

            return Fit(dist_base, samples.Select(v => (ddouble)v).ToArray(), fitting_quantile_range);
        }

        public static Distribution Fit(
            Distribution dist_base,
            ReadOnlyCollection<double> samples,
            (double min, double max) fitting_quantile_range) {

            return Fit(dist_base, samples.Select(v => (ddouble)v).ToArray(), fitting_quantile_range);
        }

        public static Distribution Fit(
            Distribution dist_base,
            IEnumerable<ddouble> samples,
            (ddouble min, ddouble max) fitting_quantile_range) {

            return FitForSortedSamples(dist_base, new ReadOnlyCollection<ddouble>(samples.Sort().ToArray()), fitting_quantile_range);
        }

        public static Distribution Fit(
            Distribution dist_base,
            ReadOnlyCollection<ddouble> samples,
            (ddouble min, ddouble max) fitting_quantile_range) {

            return FitForSortedSamples(dist_base, new ReadOnlyCollection<ddouble>(samples.Sort().ToArray()), fitting_quantile_range);
        }

        internal static Distribution FitForSortedSamples(
            Distribution dist_base,
            ReadOnlyCollection<ddouble> sorted_samples,
            (ddouble min, ddouble max) fitting_quantile_range) {
            if (!(fitting_quantile_range.min < fitting_quantile_range.max && fitting_quantile_range.min >= 0d && fitting_quantile_range.max <= 1d)) {

                throw new ArgumentException("Invalid range: min < max", nameof(fitting_quantile_range));
            }

            ddouble q_inv = 1d / (ddouble)(sorted_samples.Count - 1);

            int n = 0;
            ddouble sum_x = 0d, sum_y = 0d, sum_x2 = 0d, sum_xy = 0d;

            for (int i = 0; i < sorted_samples.Count; i++) {
                ddouble q = i * q_inv;
                if (!(q >= fitting_quantile_range.min && q <= fitting_quantile_range.max)) {
                    continue;
                }

                ddouble x = dist_base.Quantile(q), y = sorted_samples[i];

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

                return dist;
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }
    }
}
