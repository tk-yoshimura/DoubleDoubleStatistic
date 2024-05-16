using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.MultiVariateDistributions;

namespace DoubleDoubleStatistic.Misc {
    public interface IFittableMultiVariateDistribution<Distribution> where Distribution : MultiVariateDistribution {
        static abstract Distribution? Fit(IEnumerable<double[]> samples);
        static abstract Distribution? Fit(IEnumerable<ddouble[]> samples);
        static abstract Distribution? Fit(IEnumerable<Vector> samples);
    }
}
