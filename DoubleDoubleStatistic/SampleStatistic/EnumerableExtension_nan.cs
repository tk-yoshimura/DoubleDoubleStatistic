using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        public static IEnumerable<double> NaNAs(this IEnumerable<double> xs, double value) {
            return xs.Select(x => double.IsNaN(x) ? value : x);
        }

        public static IEnumerable<ddouble> NaNAs(this IEnumerable<ddouble> xs, ddouble value) {
            return xs.Select(x => ddouble.IsNaN(x) ? value : x);
        }

        public static IEnumerable<double> NaNAsZero(this IEnumerable<double> xs) {
            return xs.Select(x => double.IsNaN(x) ? 0d : x);
        }

        public static IEnumerable<ddouble> NaNAsZero(this IEnumerable<ddouble> xs) {
            return xs.Select(x => ddouble.IsNaN(x) ? 0d : x);
        }

        public static IEnumerable<double> TrimNaN(this IEnumerable<double> xs) {
            return xs.Where(x => !double.IsNaN(x));
        }

        public static IEnumerable<ddouble> TrimNaN(this IEnumerable<ddouble> xs) {
            return xs.Where(x => !ddouble.IsNaN(x));
        }
    }
}
