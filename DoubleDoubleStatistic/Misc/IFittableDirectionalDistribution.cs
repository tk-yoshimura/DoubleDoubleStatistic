using DoubleDoubleStatistic.DirectionalDistributions;

namespace DoubleDoubleStatistic.Misc {
    public interface IFittableDirectionalDistribution<Distribution, DD, D> where Distribution : DirectionalDistribution<DD, D> {
        static abstract Distribution? Fit(IEnumerable<D> samples);
        static abstract Distribution? Fit(IEnumerable<DD> samples);
    }
}
