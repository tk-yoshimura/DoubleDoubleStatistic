using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableUtil {
        public static IEnumerable<double> Linspace(double min, double max, int count, bool end_point) {
            ArgumentOutOfRangeException.ThrowIfLessThan(count, 1, nameof(count));

            double range = max - min;
            int n = end_point ? (count - 1) : count;

            for (int i = 0; i < count; i++) {
                yield return range * i / n + min;
            }
        }

        public static IEnumerable<ddouble> Linspace(ddouble min, ddouble max, int count, bool end_point) {
            ArgumentOutOfRangeException.ThrowIfLessThan(count, 1, nameof(count));

            ddouble range = max - min;
            int n = end_point ? (count - 1) : count;

            for (int i = 0; i < count; i++) {
                yield return range * i / n + min;
            }
        }
    }
}
