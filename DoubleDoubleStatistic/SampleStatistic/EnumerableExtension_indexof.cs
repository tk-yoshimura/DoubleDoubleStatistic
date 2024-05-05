using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        public static int MaxIndexOf(this IEnumerable<double> xs) {
            return xs.Select((x, index) => (x, index)).MaxBy(item => item.x).index;
        }

        public static int MaxIndexOf(this IEnumerable<ddouble> xs) {
            return xs.Select((x, index) => (x, index)).MaxBy(item => item.x).index;
        }

        public static int MinIndexOf(this IEnumerable<double> xs) {
            return xs.Select((x, index) => (x, index)).MinBy(item => item.x).index;
        }

        public static int MinIndexOf(this IEnumerable<ddouble> xs) {
            return xs.Select((x, index) => (x, index)).MinBy(item => item.x).index;
        }
    }
}
