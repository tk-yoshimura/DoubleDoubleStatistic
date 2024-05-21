using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BenfordDistribution : DiscreteDistribution {

        public int N { get; }

        private readonly ddouble logn_inv;

        private readonly Roulette roulette;

        public BenfordDistribution(int n = 10) {
            ParamAssert.ValidateShape(nameof(n), n >= 2d && n <= 1024d);

            N = n;

            logn_inv = 1d / Log(n);

            double[] probs = new double[n - 1];
            for (int k = 1; k < n; k++) {
                probs[k - 1] = (double)PMF(k);
            }
            roulette = new(probs);
        }

        public override ddouble PMF(int k) {
            return k >= 1 && k < N ? Log(1d + 1d / (ddouble)k) * logn_inv : 0d;
        }

        public override int Sample(Random random) {
            int n = roulette.NextIndex(random) + 1;

            return n;
        }

        public override (int min, int max) Support => (1, N - 1);

        public override ddouble Mean => N - LogGamma(N + 1) * logn_inv;

        public override ddouble Variance => N * N + (2d * LogBarnesG(N + 1) - (2*N - 1) * LogGamma(N + 1)) * logn_inv - Square(Mean);

        public override ddouble Skewness {
            get {
                if (N < 3) {
                    return NaN;
                }

                ddouble mu = Mean, var = Variance;

                ddouble skewness = (new int[N - 1]).Select((_, k) => PMF(k + 1) * Cube(k + 1 - mu)).Sum() / ExMath.Pow3d2(var);

                return skewness;
            }
        }
    
        public override ddouble Kurtosis {
            get {
                if (N < 4) {
                    return NaN;
                }

                ddouble mu = Mean, var = Variance;

                ddouble kurtosis = (new int[N - 1]).Select((_, k) => PMF(k + 1) * Pow(k + 1 - mu, 4)).Sum() / Square(var) - 3d;

                return kurtosis;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 1, N - 1);

        public override string Formula => "f(k; n) := log(1 + 1 / k) / log(n)";

        public override string ToString() {
            return $"{typeof(BenfordDistribution).Name}[n={N}]";
        }

    }
}
