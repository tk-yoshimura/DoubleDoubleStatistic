using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class YuleSimonDistribution : DiscreteDistribution,
        IFittableDiscreteDistribution<YuleSimonDistribution> {

        public ddouble Rho { get; }

        public YuleSimonDistribution(ddouble rho) {
            ParamAssert.ValidateShape(nameof(rho), ParamAssert.IsFinitePositive(rho));

            Rho = rho;
        }

        public override ddouble PMF(int k) {
            return k >= 1 ? Rho * Beta(k, Rho + 1d) : 0d;
        }

        public override int Sample(Random random) {
            double r1 = double.Log(random.NextUniformOpenInterval1());
            double r2 = double.Log(random.NextUniformOpenInterval1());

            double u = r1 / double.Log(1d - double.Exp(r2 / (double)Rho));

            int n = u < int.MaxValue ? int.Max(1, (int)double.Ceiling(u)) : int.MaxValue;

            return n;
        }

        public override (int min, int max) Support => (1, int.MaxValue);

        public override ddouble Mean => Rho > 1d
            ? Rho / (Rho - 1d)
            : NaN;

        public override ddouble Variance => Rho > 2d
            ? Rho * Rho / (Square(Rho - 1d) * (Rho - 2d))
            : NaN;

        public override ddouble Skewness => Rho > 3d
            ? Square(Rho + 1d) * Sqrt(Rho - 2d) / (Rho * (Rho - 3d))
            : NaN;

        public override ddouble Kurtosis => Rho > 4d
            ? (-22d + Rho * (-49d + Rho * Rho * 11d)) / (Rho * (Rho - 3d) * (Rho - 4d)) + Rho + 3d
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 1, 65536);

        public static YuleSimonDistribution? Fit(IEnumerable<int> samples) {
            if (!samples.Any() || samples.Any(n => n < 1)) {
                return null;
            }

            ddouble mean = samples.Select(n => (ddouble)n).Mean();
            ddouble pho = mean / (mean - 1d);

            try {
                return new YuleSimonDistribution(pho);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string Formula => "f(k; rho) := rho * beta(k, rho + 1)";

        public override string ToString() {
            return $"{typeof(YuleSimonDistribution).Name}[rho={Rho}]";
        }

    }
}
