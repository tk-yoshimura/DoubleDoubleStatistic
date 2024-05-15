using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class PoissonDistribution : DiscreteDistribution,
        IFittableDiscreteDistribution<PoissonDistribution> {

        const int random_gen_max_index = 65536;

        public ddouble Lambda { get; }

        private readonly ddouble lnlambda;

        private Roulette? roulette = null;

        public PoissonDistribution(ddouble lambda) {
            ValidateScale(lambda);

            Lambda = lambda;
            lnlambda = Log(lambda);
        }

        public override ddouble PMF(int k) {
            return k >= 0 ? Exp(-LogGamma(k + 1) + lnlambda * k - Lambda) : 0d;
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

        public override ddouble Mean => Lambda;

        public override ddouble Variance => Lambda;

        public override ddouble Skewness => 1d / Sqrt(Lambda);

        public override ddouble Kurtosis => 1d / Lambda;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, 65535);

        public static PoissonDistribution? Fit(IEnumerable<int> samples) {
            if (samples.Count() < 1 || samples.Any(n => n < 0)) {
                return null;
            }

            ddouble mean = samples.Select(n => (ddouble)n).Mean();

            try {
                return new PoissonDistribution(mean);
            }
            catch {
                return null;
            }
        }

        public override string Formula => "f(k; lambda) := lambda^k * exp(-lambda) / k!";

        public override string ToString() {
            return $"{typeof(PoissonDistribution).Name}[lambda={Lambda}]";
        }

    }
}
