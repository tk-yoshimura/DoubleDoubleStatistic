using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BinomialDistribution : DiscreteDistribution {

        public int N { get; }
        public ddouble P { get; }

        private readonly ddouble p, q;

        private const int random_gen_batches = 32;
        private readonly Roulette? roulette_x32 = null, roulette_r32 = null;

        public BinomialDistribution(int n, ddouble p) {
            ValidateShape(n, n => n > 0);
            ValidateShape(p, p => p >= 0d && p <= 1d);

            N = n;
            P = p;
            this.p = p;
            q = 1d - p;

            if (n >= random_gen_batches) {
                roulette_x32 = new Roulette(new double[random_gen_batches + 1].Select(
                    (_, k) => (double)(Binomial(random_gen_batches, k) * Pow(p, k) * Pow(q, random_gen_batches - k)
                )));
            }
            if ((n % random_gen_batches) > 0) {
                int r = n % random_gen_batches;

                roulette_r32 = new Roulette(new double[r + 1].Select(
                    (_, k) => (double)(Binomial(r, k) * Pow(p, k) * Pow(q, r - k)
                )));
            }
        }

        public override ddouble PMF(int k) {
            return (k >= 0 && k <= N) ? Binomial(N, k) * Pow(p, k) * Pow(q, N - k) : 0d;
        }

        public override int Sample(Random random) {
            int i = N, cnt = 0;

            if (roulette_x32 is not null) {
                while (i >= random_gen_batches) {
                    cnt += roulette_x32.NextIndex(random);
                    i -= random_gen_batches;
                }
            }

            if (i > 0 && roulette_r32 is not null) {
                cnt += roulette_r32.NextIndex(random);
            }

            return cnt;
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
