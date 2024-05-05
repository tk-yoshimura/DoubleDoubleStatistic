using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LogarithmicDistribution : DiscreteDistribution {

        public ddouble P { get; }

        private readonly ddouble p, q, lnq;

        public LogarithmicDistribution(ddouble p) {
            ValidateShape(p, p => p > 0d && p < 1d);

            P = p;

            this.p = p;
            q = 1d - p;
            lnq = Log(q);
        }

        public override ddouble PMF(long k) {
            return k >= 1 ? -Pow(p, k) / (lnq * k) : 0d;
        }

        public override (long min, long max) Support => (1, long.MaxValue);

        public override ddouble Mean => -p / (lnq * q);

        public override ddouble Variance => -p * (p + lnq) / Square(q * lnq);

        public override ddouble Skewness => -lnq * (p * (3d * lnq + 2d * p) + (1d + p) * lnq * lnq) / (lnq * (p + lnq) * Sqrt(-p * (p + lnq)));

        public override ddouble Kurtosis => -(6d * p * p * p + 12d * p * p * lnq + p * (4d * p + 7d) * lnq * lnq + (p * p + 4d * p + 1d) * lnq * lnq * lnq) / (p * Square(p + lnq));

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 1, 65536);

        public override string Formula => "f(k; p) := -1 / log(1 - p) * p^k / k";

        public override string ToString() {
            return $"{typeof(LogarithmicDistribution).Name}[p={P}]";
        }

    }
}
