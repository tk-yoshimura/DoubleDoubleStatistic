using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic.InternalUtils {
    public static class Indexer {
        public static int BisectionSearch(double x, ReadOnlyCollection<double> thresholds) {
            Debug.Assert(thresholds.Count >= 2);

            if (x >= thresholds[^1]) {
                return thresholds.Count - 2;
            }

            int index = 0;

            for (int h = thresholds.Count / 2; h >= 1; h /= 2) {
                for (int i = index, m = thresholds.Count - h; i < m; i += h) {
                    if (x < thresholds[i + h]) {
                        index = i;
                        break;
                    }
                }
            }

            return index;
        }

        public static int BisectionSearch(ddouble x, ReadOnlyCollection<ddouble> thresholds) {
            Debug.Assert(thresholds.Count >= 2);

            if (x >= thresholds[^1]) {
                return thresholds.Count - 2;
            }

            int index = 0;

            for (int h = thresholds.Count / 2; h >= 1; h /= 2) {
                for (int i = index, m = thresholds.Count - h; i < m; i += h) {
                    if (x < thresholds[i + h]) {
                        index = i;
                        break;
                    }
                }
            }

            return index;
        }
    }
}
