using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class YuleSimonDistribution : DiscreteDistribution {

        public ddouble Rho { get; }

        public YuleSimonDistribution(ddouble rho) {
            ValidateShape(rho, rho => rho > 0d);

            Rho = rho;
        }

        public override ddouble PMF(long k) {
            return k >= 1 ? Rho * Beta(k, Rho + 1d) : 0d;
        }

        public override (long min, long max) Support => (1, long.MaxValue);

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

        public override string Formula => "f(k; rho) := rho * beta(k, rho + 1)";

        public override string ToString() {
            return $"{typeof(YuleSimonDistribution).Name}[rho={Rho}]";
        }

    }
}
