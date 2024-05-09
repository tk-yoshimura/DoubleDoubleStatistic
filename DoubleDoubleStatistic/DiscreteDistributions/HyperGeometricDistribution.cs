using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class HyperGeometricDistribution : DiscreteDistribution {

        public int N { get; }
        public int M { get; }
        public int R { get; }

        private readonly ddouble binom_nr_inv;
        private readonly int min_k;
        private readonly Roulette roulette;

        public HyperGeometricDistribution(int n, int m, int r) {
            ValidateShape(n, n => n > 0);
            ValidateShape(m, m => m >= 0 && m <= n);
            ValidateShape(r, r => r >= 0 && r <= n);

            N = n;
            M = m;
            R = r;

            binom_nr_inv = 1d / Binomial(n, r);
            min_k = Support.min;

            double[] probs = new double[Support.max - Support.min + 1];
            for (int i = 0; i < probs.Length; i++) {
                probs[i] = (double)PMF(i + min_k);
            }

            roulette = new(probs);
        }

        public override ddouble PMF(int k) {
            return k >= Support.min && k <= Support.max ? Binomial(M, k) * Binomial(N - M, R - k) * binom_nr_inv : 0d;
        }

        public override int Sample(Random random) {
            return roulette.NextIndex(random) + min_k;
        }

        public override (int min, int max) Support => (int.Max(0, checked(R + M - N)), int.Min(R, M));

        public override ddouble Mean => R * M / (ddouble)N;

        public override ddouble Variance => N > 1
            ? (ddouble)checked(R * M * (N - M) * (N - R)) / checked(N * N * (N - 1))
            : NaN;

        public override ddouble Skewness => N > 2
            ? checked((N - 2 * M) * (N - 2 * R)) * Sqrt(N - 1) / ((N - 2) * Sqrt(checked(R * M * (N - M) * (N - R))))
            : 0;

        public override ddouble Kurtosis => N > 3
            ? (ddouble)checked((N - 1) * N * N * (N * (N + 1) - 6 * M * (N - M) - 6 * R * (N - R)) + 6 * R * M * (N - M) * (N - R) * (5 * N - 6))
              / checked(R * M * (N - M) * (N - R) * (N - 2) * (N - 3))
            : NaN;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, Support.min, Support.max);

        public override string Formula => "f(k; n, m, r) := binom(m, k) * binom(n - m, r - k) / binom(n, r)";

        public override string ToString() {
            return $"{typeof(HyperGeometricDistribution).Name}[n={N},m={M},r={R}]";
        }

    }
}
