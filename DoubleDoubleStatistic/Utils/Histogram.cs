using DoubleDouble;
using System.Collections.ObjectModel;

namespace DoubleDoubleStatistic.Utils {
    public class Histogram {

        public int Bins { get; }

        public ReadOnlyCollection<ddouble> BinsEdge { get; }
        public ReadOnlyCollection<ddouble> BinsCentroid { get; }

        public ReadOnlyCollection<int> Counts { get; }
        public ReadOnlyCollection<ddouble> Density { get; }

        public Histogram(IEnumerable<double> samples, int bins, (double min, double max) range)
            : this(new ReadOnlyCollection<double>(samples.ToArray()), bins, range) { }

        public Histogram(IEnumerable<ddouble> samples, int bins, (ddouble min, ddouble max) range)
            : this(new ReadOnlyCollection<ddouble>(samples.ToArray()), bins, range) { }

        public Histogram(ReadOnlyCollection<double> samples, int bins, (double min, double max) range)
            : this(new ReadOnlyCollection<ddouble>(samples.Select(v => (ddouble)v).ToArray()), bins, ((ddouble)range.min, (ddouble)range.max)) { }

        public Histogram(ReadOnlyCollection<ddouble> samples, int bins, (ddouble min, ddouble max) range) {
            ArgumentOutOfRangeException.ThrowIfLessThan(bins, 1, nameof(bins));
            if (!(ddouble.IsFinite(range.min) && ddouble.IsFinite(range.max) && range.min < range.max)) {
                throw new ArgumentException("invalid range: min < max", nameof(range));
            }

            ddouble xr = range.max - range.min, h = xr / bins;

            int[] counts = new int[bins];

            foreach (ddouble sample in samples) {
                if (!(sample >= range.min && sample <= range.max)) {
                    continue;
                }

                int index = int.Min(bins - 1, (int)ddouble.Clamp(ddouble.Floor((sample - range.min) / h), 0d, bins));

                counts[index]++;
            }

            ddouble w = 1d / (h * samples.Count);

            ddouble[] density = counts.Select(n => n * w).ToArray();

            Bins = bins;
            BinsEdge = new ReadOnlyCollection<ddouble>((new ddouble[bins + 1]).Select((_, idx) => range.min + xr * idx / bins).ToArray());
            BinsCentroid = new ReadOnlyCollection<ddouble>((new ddouble[bins]).Select((_, idx) => range.min + xr * (idx + 0.5d) / bins).ToArray());
            Counts = new ReadOnlyCollection<int>(counts);
            Density = new ReadOnlyCollection<ddouble>(density);
        }
    }
}
