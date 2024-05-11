using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System;
using System.Collections.ObjectModel;

namespace DoubleDoubleStatistic.Utils {
    public class Histogram {

        public int Bins { get; }

        public ReadOnlyCollection<ddouble> BinsEdge { get; }
        public ReadOnlyCollection<ddouble> BinsCentor { get; }

        public ReadOnlyCollection<int> Counts { get; }
        public ReadOnlyCollection<ddouble> Density { get; }

        public Histogram(IEnumerable<double> samples, int bins, (double min, double max) range)
            : this(new ReadOnlyCollection<double>(samples.ToArray()), bins, range) { }

        public Histogram(ReadOnlyCollection<double> samples, int bins, (double min, double max) range)
            : this(new ReadOnlyCollection<ddouble>(samples.Select(v => (ddouble)v).ToArray()), bins, ((ddouble)range.min, (ddouble)range.max)) { }

        public Histogram(IEnumerable<ddouble> samples, int bins, (ddouble min, ddouble max) range)
            : this(new ReadOnlyCollection<ddouble>(samples.ToArray()), bins, range) { }

        public Histogram(ReadOnlyCollection<ddouble> samples, int bins, (ddouble min, ddouble max) range) {
            ArgumentOutOfRangeException.ThrowIfLessThan(bins, 1, nameof(bins));
            if (!(ddouble.IsFinite(range.min) && ddouble.IsFinite(range.max) && range.min < range.max)) {
                throw new ArgumentException("Invalid range: min < max", nameof(range));
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
            BinsCentor = new ReadOnlyCollection<ddouble>((new ddouble[bins]).Select((_, idx) => range.min + xr * (idx + 0.5d) / bins).ToArray());
            Counts = new ReadOnlyCollection<int>(counts);
            Density = new ReadOnlyCollection<ddouble>(density);
        }

        public Histogram(IEnumerable<double> samples, IEnumerable<double> bins_edge)
            : this(
                  new ReadOnlyCollection<ddouble>(samples.Select(v => (ddouble)v).ToArray()), 
                  new ReadOnlyCollection<ddouble>(bins_edge.Select(v => (ddouble)v).ToArray())
            ) {}

        public Histogram(IEnumerable<double> samples, ReadOnlyCollection<double> bins_edge)
            : this(samples, (IEnumerable<double>)bins_edge) {}

        public Histogram(ReadOnlyCollection<double> samples, IEnumerable<double> bins_edge)
            : this((IEnumerable<double>)samples, bins_edge) {}

        public Histogram(ReadOnlyCollection<double> samples, ReadOnlyCollection<double> bins_edge)
            : this((IEnumerable<double>)samples, (IEnumerable<double>)bins_edge) {}

        public Histogram(IEnumerable<ddouble> samples, IEnumerable<ddouble> bins_edge)
            : this(new ReadOnlyCollection<ddouble>(samples.ToArray()), new ReadOnlyCollection<ddouble>(bins_edge.ToArray())) {}

        public Histogram(IEnumerable<ddouble> samples, ReadOnlyCollection<ddouble> bins_edge)
            : this(new ReadOnlyCollection<ddouble>(samples.ToArray()), bins_edge) {}

        public Histogram(ReadOnlyCollection<ddouble> samples, IEnumerable<ddouble> bins_edge)
            : this(samples, new ReadOnlyCollection<ddouble>(bins_edge.ToArray())) {}

        public Histogram(ReadOnlyCollection<ddouble> samples, ReadOnlyCollection<ddouble> bins_edge) {
            ArgumentOutOfRangeException.ThrowIfLessThan(bins_edge.Count, 2, nameof(bins_edge));
            for (int i = 0; i < bins_edge.Count - 1; i++) {
                if (!(bins_edge[i] < bins_edge[i + 1])) {
                    throw new ArgumentException("Invalid bins: sequential", nameof(bins_edge));
                }
            }
            if (!bins_edge.All(ddouble.IsFinite)) { 
                throw new ArgumentException("Invalid bins: finite", nameof(bins_edge));
            }

            ddouble xmin = bins_edge[0], xmax = bins_edge[^1];

            int[] counts = new int[bins_edge.Count - 1];

            foreach (ddouble sample in samples) {
                if (!(sample >= xmin && sample <= xmax)) {
                    continue;
                }

                int index = Indexer.BisectionSearch(sample, bins_edge);

                counts[index]++;
            }

            ddouble[] density = counts.Select((n, idx) => n / ((bins_edge[idx + 1] - bins_edge[idx]) * samples.Count)).ToArray();

            Bins = bins_edge.Count - 1;
            BinsEdge = new ReadOnlyCollection<ddouble>(bins_edge);
            BinsCentor = new ReadOnlyCollection<ddouble>(counts.Select((_, idx) => (bins_edge[idx + 1] + bins_edge[idx]) * 0.5d).ToArray());
            Counts = new ReadOnlyCollection<int>(counts);
            Density = new ReadOnlyCollection<ddouble>(density);
        }
    }
}
