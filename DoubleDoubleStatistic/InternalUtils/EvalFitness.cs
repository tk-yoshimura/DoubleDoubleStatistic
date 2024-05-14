using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;

namespace DoubleDoubleStatistic.InternalUtils {
    internal static class EvalFitness {
        public static ddouble MeanSquaredError(ContinuousDistribution dist, ddouble[] qs, ddouble[] ys) {
            Debug.Assert(qs.Length == ys.Length);

            ddouble error = qs.Select(q => dist.Quantile(q)).Select((x, idx) => ddouble.Square(x - ys[idx])).Mean();

            return error;
        }
    }
}
