using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;

namespace DoubleDoubleStatistic.Misc {
    public interface IFittableContinuousDistribution<Distribution> where Distribution : ContinuousDistribution {
        static abstract (Distribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions);

        static abstract (Distribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions);
    }
}
