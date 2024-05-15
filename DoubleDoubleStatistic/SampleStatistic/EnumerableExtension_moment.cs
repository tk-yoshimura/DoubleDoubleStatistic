using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {
        public static double Mean(this IEnumerable<double> xs) {
            if (!xs.Any()) {
                return double.NaN;
            }

            return xs.Average();
        }

        public static ddouble Mean(this IEnumerable<ddouble> xs) {
            if (!xs.Any()) {
                return double.NaN;
            }

            return xs.Average();
        }

        public static double Variance(this IEnumerable<double> xs) {
            double mean = xs.Mean();
            double variance = xs.Select(x => { double d = x - mean; return d * d; }).Mean();

            return variance;
        }

        public static ddouble Variance(this IEnumerable<ddouble> xs) {
            ddouble mean = xs.Mean();
            ddouble variance = xs.Select(x => { ddouble d = x - mean; return d * d; }).Mean();

            return variance;
        }

        public static double Skewness(this IEnumerable<double> xs) {
            (double mean, double variance) = xs.MeanVariance();
            double skewness = xs.Select(x => { double d = x - mean; return d * d * d; }).Mean() / double.Pow(double.Sqrt(variance), 3);

            return skewness;
        }

        public static ddouble Skewness(this IEnumerable<ddouble> xs) {
            (ddouble mean, ddouble variance) = xs.MeanVariance();
            ddouble skewness = xs.Select(x => { ddouble d = x - mean; return d * d * d; }).Mean() / ExMath.Pow3d2(variance);

            return skewness;
        }

        public static double Kurtosis(this IEnumerable<double> xs) {
            (double mean, double variance) = xs.MeanVariance();
            double kurtosis = xs.Select(x => { double d = x - mean; return d * d * d * d; }).Mean() / (variance * variance) - 3d;

            return kurtosis;
        }

        public static ddouble Kurtosis(this IEnumerable<ddouble> xs) {
            (ddouble mean, ddouble variance) = xs.MeanVariance();
            ddouble kurtosis = xs.Select(x => { ddouble d = x - mean; return d * d * d * d; }).Mean() / (variance * variance) - 3d;

            return kurtosis;
        }

        public static (double mean, double variance) MeanVariance(this IEnumerable<double> xs) {
            double mean = xs.Mean();
            double variance = xs.Select(x => { double d = x - mean; return d * d; }).Mean();

            return (mean, variance);
        }

        public static (ddouble mean, ddouble variance) MeanVariance(this IEnumerable<ddouble> xs) {
            ddouble mean = xs.Mean();
            ddouble variance = xs.Select(x => { ddouble d = x - mean; return d * d; }).Mean();

            return (mean, variance);
        }

        public static double StandardDeviation(this IEnumerable<double> xs) {
            return double.Sqrt(xs.Variance());
        }

        public static ddouble StandardDeviation(this IEnumerable<ddouble> xs) {
            return ddouble.Sqrt(xs.Variance());
        }
    }
}
