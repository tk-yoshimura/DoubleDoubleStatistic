using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        public static IEnumerable<double> Sort(this IEnumerable<double> xs) {
            return xs.OrderBy(x => x);
        }

        public static IEnumerable<ddouble> Sort(this IEnumerable<ddouble> xs) {
            return xs.OrderBy(x => x);
        }
    }
}
