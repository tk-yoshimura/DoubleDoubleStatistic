using DoubleDouble;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        public static IEnumerable<double> Trim(this IEnumerable<double> xs, double quantile_min, double quantile_max) {
            if (quantile_min > 0d && quantile_max < 1d) {
                (double xmin, double xmax) = xs.Quantile(quantile_min, quantile_max);

                return xs.Where(x => x >= xmin && x <= xmax);
            }

            if (quantile_min > 0d) {
                double xmin = xs.Quantile(quantile_min);

                return xs.Where(x => x >= xmin);
            }

            if (quantile_max < 1d) {
                double xmax = xs.Quantile(quantile_max);

                return xs.Where(x => x <= xmax);
            }

            return xs;
        }

        public static IEnumerable<ddouble> Trim(this IEnumerable<ddouble> xs, ddouble quantile_min, ddouble quantile_max) {
            if (quantile_min > 0d && quantile_max < 1d) {
                (ddouble xmin, ddouble xmax) = xs.Quantile(quantile_min, quantile_max);

                return xs.Where(x => x >= xmin && x <= xmax);
            }

            if (quantile_min > 0d) {
                ddouble xmin = xs.Quantile(quantile_min);

                return xs.Where(x => x >= xmin);
            }

            if (quantile_max < 1d) {
                ddouble xmax = xs.Quantile(quantile_max);

                return xs.Where(x => x <= xmax);
            }

            return xs;
        }
    }
}
