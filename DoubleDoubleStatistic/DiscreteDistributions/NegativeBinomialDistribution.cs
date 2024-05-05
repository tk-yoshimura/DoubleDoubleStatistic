using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class NegativeBinomialDistribution : DiscreteDistribution {

        public int N { get; }
        public ddouble P { get; }

        private readonly ddouble p, q, p_pown;

        public NegativeBinomialDistribution(int n, ddouble p) {
            ValidateShape(n, n => n > 0);
            ValidateShape(p, p => p >= 0d && p <= 1d);

            N = n;
            P = p;
            this.p = p;
            q = 1d - p;
            p_pown = Pow(p, N);
        }

        public override ddouble PMF(long k) {
            return k >= 0 ? Binomial((int)k + N - 1, (int)k) * p_pown * Pow(q, k) : 0d;
        }

        public override ddouble Mean => N * q / p;

        public override ddouble Variance => N * q / (p * p);

        public override ddouble Skewness => (2d - p) / Sqrt(N * q);

        public override ddouble Kurtosis => (p * p) / (N * q) + 6d / N;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, 65535);

        public override string Formula => "f(k; n, p) := binom(k + n - 1, k) * p^N * q^k";

        public override string ToString() {
            return $"{typeof(NegativeBinomialDistribution).Name}[n={N},p={P}]";
        }

    }
}
