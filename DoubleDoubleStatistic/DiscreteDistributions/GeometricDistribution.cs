using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class GeometricDistribution : DiscreteDistribution {

        public ddouble P { get; }

        private readonly ddouble p, q;

        public GeometricDistribution(ddouble p) {
            ValidateShape(p, p => p > 0d && p <= 1d);

            P = p;

            this.p = p;
            q = 1d - p;
        }

        public override ddouble PMF(long k) {
            return k >= 0 ? p * Pow(q, k) : 0d;
        }

        public override ddouble Mean => q / p;

        public override ddouble Variance => q / Square(p);

        public override ddouble Skewness => (2d - p) / Sqrt(q);

        public override ddouble Kurtosis => Square(p) / q + 6d;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, 65535);

        public override string Formula => "f(k; p) := p * (1 - p)^k";

        public override string ToString() {
            return $"{typeof(GeometricDistribution).Name}[p={P}]";
        }

    }
}
