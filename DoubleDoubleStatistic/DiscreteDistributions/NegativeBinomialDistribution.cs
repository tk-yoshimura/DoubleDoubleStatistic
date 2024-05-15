using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class NegativeBinomialDistribution : DiscreteDistribution,
        IFittableDiscreteDistribution<NegativeBinomialDistribution> {

        const int random_gen_max_index = 65536;

        public int N { get; }
        public ddouble P { get; }

        private readonly ddouble p, q, p_pown;

        private Roulette? roulette = null;

        public NegativeBinomialDistribution(int n, ddouble p) {
            ValidateShape(n, n => n > 0);
            ValidateShape(p, p => p >= 0d && p <= 1d);

            N = n;
            P = p;
            this.p = p;
            q = 1d - p;
            p_pown = Pow(p, N);
        }

        public override ddouble PMF(int k) {
            return k >= 0 ? Binomial(k + N - 1, k) * p_pown * Pow(q, k) : 0d;
        }

        public override int Sample(Random random) {
            roulette ??= SetupRoulette();

            return roulette.NextIndex(random);
        }

        public Roulette SetupRoulette() {
            List<double> probs = [];
            ddouble ps = 0d;
            for (int k = 0; k <= random_gen_max_index && (double)ps < 1d; k++) {
                ddouble p = PMF(k);
                ps += p;

                probs.Add((double)p);
            }

            return new(probs);
        }

        public override ddouble Mean => N * q / p;

        public override ddouble Variance => N * q / (p * p);

        public override ddouble Skewness => (2d - p) / Sqrt(N * q);

        public override ddouble Kurtosis => (p * p) / (N * q) + 6d / N;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, 65535);

        public static NegativeBinomialDistribution? Fit(IEnumerable<int> samples) {
            if (!samples.Any() || samples.Any(n => n < 0)) {
                return null;
            }

            (ddouble mean, ddouble variance) = samples.Select(n => (ddouble)n).MeanVariance();
            int N = (int)Clamp(Round(mean * mean / (variance - mean)), 1, int.MaxValue);
            ddouble p = N / (N + mean);

            try {
                return new NegativeBinomialDistribution(N, p);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string Formula => "f(k; n, p) := binom(k + n - 1, k) * p^N * q^k";

        public override string ToString() {
            return $"{typeof(NegativeBinomialDistribution).Name}[n={N},p={P}]";
        }

    }
}
