using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        private static double QuantileIndexer(List<double> sorted_list, double q) {
            int count = sorted_list.Count;

            if (count < 1 || !(q >= 0d && q <= 1d)) {
                return double.NaN;
            }

            if (q == 0d) {
                return sorted_list[0];
            }

            if (q == 1d) {
                return sorted_list[^1];
            }

            double index = double.Clamp((count - 1) * q, 0, count - 1);

            int index_0 = int.Max(0, (int)double.Floor(index)), index_1 = int.Min(count - 1, (int)double.Floor(index + 1));
            double c = index - index_0;

            double x = sorted_list[index_0] * (1d - c) + sorted_list[index_1] * c;

            return x;
        }

        private static ddouble QuantileIndexer(List<ddouble> sorted_list, ddouble q) {
            int count = sorted_list.Count;

            if (count < 1 || !(q >= 0d && q <= 1d)) {
                return ddouble.NaN;
            }

            if (q == 0d) {
                return sorted_list[0];
            }

            if (q == 1d) {
                return sorted_list[^1];
            }

            ddouble index = ddouble.Clamp((count - 1) * q, 0, count - 1);

            int index_0 = int.Max(0, (int)ddouble.Floor(index)), index_1 = int.Min(count - 1, (int)ddouble.Floor(index + 1));
            ddouble c = index - index_0;

            ddouble x = sorted_list[index_0] * (1d - c) + sorted_list[index_1] * c;

            return x;
        }

        public static double Quantile(this IEnumerable<double> xs, double q) {
            List<double> sorted_list = [.. xs.Sort()];

            return QuantileIndexer(sorted_list, q);
        }

        public static ddouble Quantile(this IEnumerable<ddouble> xs, ddouble q) {
            List<ddouble> sorted_list = [.. xs.Sort()];

            return QuantileIndexer(sorted_list, q);
        }

        public static IEnumerable<double> Quantile(this IEnumerable<double> xs, params double[] qs) {
            List<double> sorted_list = [.. xs.Sort()];

            foreach (double q in qs) {
                yield return QuantileIndexer(sorted_list, q);
            }
        }

        public static IEnumerable<ddouble> Quantile(this IEnumerable<ddouble> xs, params ddouble[] qs) {
            List<ddouble> sorted_list = [.. xs.Sort()];

            foreach (ddouble q in qs) {
                yield return QuantileIndexer(sorted_list, q);
            }
        }

        public static (double x0, double x1) Quantile(this IEnumerable<double> xs, double q0, double q1) {
            double[] ys = xs.Quantile([q0, q1]).ToArray();

            return (ys[0], ys[1]);
        }

        public static (ddouble x0, ddouble x1) Quantile(this IEnumerable<ddouble> xs, ddouble q0, ddouble q1) {
            ddouble[] ys = xs.Quantile([q0, q1]).ToArray();

            return (ys[0], ys[1]);
        }

        public static (double x0, double x1, double x2) Quantile(this IEnumerable<double> xs, double q0, double q1, double q2) {
            double[] ys = xs.Quantile([q0, q1, q2]).ToArray();

            return (ys[0], ys[1], ys[2]);
        }

        public static (ddouble x0, ddouble x1, ddouble x2) Quantile(this IEnumerable<ddouble> xs, ddouble q0, ddouble q1, ddouble q2) {
            ddouble[] ys = xs.Quantile([q0, q1, q2]).ToArray();

            return (ys[0], ys[1], ys[2]);
        }

        public static (double x0, double x1, double x2, double x3) Quantile(this IEnumerable<double> xs, double q0, double q1, double q2, double q3) {
            double[] ys = xs.Quantile([q0, q1, q2, q3]).ToArray();

            return (ys[0], ys[1], ys[2], ys[3]);
        }

        public static (ddouble x0, ddouble x1, ddouble x2, ddouble x3) Quantile(this IEnumerable<ddouble> xs, ddouble q0, ddouble q1, ddouble q2, ddouble q3) {
            ddouble[] ys = xs.Quantile([q0, q1, q2, q3]).ToArray();

            return (ys[0], ys[1], ys[2], ys[3]);
        }
    }
}
