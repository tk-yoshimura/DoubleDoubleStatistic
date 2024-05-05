using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class PoissonDistribution : DiscreteDistribution {

        public ddouble Lambda { get; }

        private readonly ddouble lnlambda;

        public PoissonDistribution(ddouble lambda) {
            ValidateScale(lambda);

            Lambda = lambda;
            lnlambda = Log(lambda);
        }

        public override ddouble PMF(long k) {
            return k >= 0 ? Exp(-LogGamma(k + 1) + lnlambda * k - Lambda) : 0d;
        }

        public override ddouble Mean => Lambda;

        public override ddouble Variance => Lambda;

        public override ddouble Skewness => 1d / Sqrt(Lambda);

        public override ddouble Kurtosis => 1d / Lambda;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, 65535);

        public override string Formula => "f(k; lambda) := lambda^k * exp(-lambda) / k!";

        public override string ToString() {
            return $"{typeof(PoissonDistribution).Name}[lambda={Lambda}]";
        }

    }
}
