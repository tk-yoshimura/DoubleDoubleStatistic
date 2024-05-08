using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BinomialDistribution : DiscreteDistribution {

        public int N { get; }
        public ddouble P { get; }

        private readonly ddouble p, q;

        public BinomialDistribution(int n, ddouble p) {
            ValidateShape(n, n => n > 0);
            ValidateShape(p, p => p >= 0d && p <= 1d);

            N = n;
            P = p;
            this.p = p;
            q = 1d - p;
        }

        public override ddouble PMF(int k) {
            return (k >= 0 && k <= N) ? Binomial(N, k) * Pow(p, k) * Pow(q, N - k) : 0d;
        }

        public override (int min, int max) Support => (0, N);

        public override ddouble Mean => N * p;

        public override ddouble Variance => N * p * q;

        public override ddouble Skewness => (q - p) / Sqrt(N * p * q);

        public override ddouble Kurtosis => 1d / (N * p * q) - 6d / N;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, N);

        public override string Formula => "f(k; n, p) := binom(n, k) * p^k * q^(n - k)";

        public override string ToString() {
            return $"{typeof(BinomialDistribution).Name}[n={N},p={P}]";
        }

    }
}
