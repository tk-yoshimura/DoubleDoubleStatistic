using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using System.Collections.ObjectModel;

namespace DoubleDoubleStatistic.Utils {
    public class QQPlot {
        public ReadOnlyCollection<ddouble> Expected { get; }
        public ReadOnlyCollection<ddouble> Sample { get; }

        public int Count => Expected.Count;

        public QQPlot(ContinuousDistribution dist, IEnumerable<double> samples)
            : this(dist, new ReadOnlyCollection<ddouble>(samples.Select(v => (ddouble)v).ToArray())) { }

        public QQPlot(ContinuousDistribution dist, ReadOnlyCollection<double> samples)
            : this(dist, new ReadOnlyCollection<ddouble>(samples.Select(v => (ddouble)v).ToArray())) { }

        public QQPlot(ContinuousDistribution dist, IEnumerable<ddouble> samples)
            : this(dist, new ReadOnlyCollection<ddouble>(samples.ToArray())) { }

        public QQPlot(ContinuousDistribution dist, ReadOnlyCollection<ddouble> samples) {
            List<ddouble> sorted_samples = samples.Sort().ToList();

            ddouble q_inv = 1d / (ddouble)(sorted_samples.Count - 1);

            List<ddouble> xs = [], ys = [];

            for (int i = 0; i < sorted_samples.Count; i++) {
                ddouble q = i * q_inv;

                ddouble x = dist.Quantile(q), y = sorted_samples[i];

                xs.Add(x);
                ys.Add(y);
            }

            Expected = new ReadOnlyCollection<ddouble>(xs);
            Sample = new ReadOnlyCollection<ddouble>(ys);
        }

        internal static ddouble MSE(
            ContinuousDistribution dist,
            ReadOnlyCollection<ddouble> sorted_samples,
            (ddouble min, ddouble max) fitting_quantile_range) {

            if (!(fitting_quantile_range.min < fitting_quantile_range.max && fitting_quantile_range.min >= 0d && fitting_quantile_range.max <= 1d)) {
                throw new ArgumentException("Invalid range: min < max", nameof(fitting_quantile_range));
            }

            ddouble q_inv = 1d / (ddouble)(sorted_samples.Count - 1);

            int n = 0;
            ddouble error = 0d;

            for (int i = 0; i < sorted_samples.Count; i++) {
                ddouble q = i * q_inv;
                if (!(q >= fitting_quantile_range.min && q <= fitting_quantile_range.max)) {
                    continue;
                }

                ddouble x = dist.Quantile(q), y = sorted_samples[i];

                n++;
                error += ddouble.Square(x - y);
            }

            error /= n;

            return error;
        }
    }
}
