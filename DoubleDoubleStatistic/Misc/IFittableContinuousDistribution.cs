using DoubleDouble;

namespace DoubleDoubleStatistic.Misc {
    public interface IFittableDistribution<Distribution> {
        static abstract (Distribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions);

        static abstract (Distribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions);
    }
}
