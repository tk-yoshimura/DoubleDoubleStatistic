using DoubleDoubleStatistic.DiscreteDistributions;

namespace DoubleDoubleStatistic.Misc {
    public interface IFittableDiscreteDistribution<Distribution> where Distribution : DiscreteDistribution {
        static abstract Distribution? Fit(IEnumerable<int> samples);
    }
}
