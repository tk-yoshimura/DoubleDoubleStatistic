using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        public static double Median(this IEnumerable<double> xs) {
            if (xs.Count() < 1) {
                return double.NaN;
            }

            List<double> sorted = [.. xs.Sort()];

            int count = sorted.Count, index = (count - 1) / 2;

            if ((count & 1) == 0) {
                return (sorted[index] + sorted[index + 1]) * 0.5d;
            }
            else {
                return sorted[index];
            }
        }

        public static ddouble Median(this IEnumerable<ddouble> xs) {
            if (xs.Count() < 1) {
                return double.NaN;
            }

            List<ddouble> sorted = [.. xs.Sort()];

            int count = sorted.Count, index = (count - 1) / 2;

            if ((count & 1) == 0) {
                return (sorted[index] + sorted[index + 1]) * 0.5d;
            }
            else {
                return sorted[index];
            }
        }
    }
}
