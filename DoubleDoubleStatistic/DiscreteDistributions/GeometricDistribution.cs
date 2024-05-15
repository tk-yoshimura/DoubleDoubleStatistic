using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class GeometricDistribution : DiscreteDistribution,
        IFittableDiscreteDistribution<GeometricDistribution> {

        public ddouble P { get; }

        private readonly ddouble p, q;

        private const int random_gen_batches = 32, random_gen_max_iter = 65536;
        private readonly Roulette? roulette_x32 = null;

        public GeometricDistribution(ddouble p) {
            ValidateShape(p, p => p > 0d && p <= 1d);

            P = p;

            this.p = p;
            q = 1d - p;

            if (p > 0.25d) {
                double[] probs = new double[random_gen_batches + 1];
                ddouble v = p, sum_v = 0;
                for (int i = 0; i < random_gen_batches; i++) {
                    probs[i] = (double)v;
                    sum_v += v;
                    v *= q;
                }

                probs[^1] = double.Max(0d, (double)(1d - sum_v));

                roulette_x32 = new(probs);
            }
        }

        public override ddouble PMF(int k) {
            return k >= 0 ? p * Pow(q, k) : 0d;
        }

        public override int Sample(Random random) {
            if (roulette_x32 is not null) {
                int cnt = 0;

                for (int i = 0; i < random_gen_max_iter; i++) {
                    int n = roulette_x32.NextIndex(random);
                    cnt += n;
                    if (n < random_gen_batches) {
                        break;
                    }
                }

                return cnt;
            }
            else {
                double u = double.Log(random.NextUniformOpenInterval0()) / double.Log((double)q);

                if (u < int.MaxValue) {
                    return (int)double.Floor(u);
                }
                else {
                    return int.MaxValue;
                }
            }
        }

        public override ddouble Mean => q / p;

        public override ddouble Variance => q / Square(p);

        public override ddouble Skewness => (2d - p) / Sqrt(q);

        public override ddouble Kurtosis => Square(p) / q + 6d;

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 0, 65535);

        public static GeometricDistribution? Fit(IEnumerable<int> samples) {
            if (samples.Count() < 1 || samples.Any(n => n < 0)) {
                return null;
            }

            ddouble mean = samples.Select(n => (ddouble)n).Mean();

            ddouble p = 1d / (mean + 1d);

            return new GeometricDistribution(p);
        }

        public override string Formula => "f(k; p) := p * (1 - p)^k";

        public override string ToString() {
            return $"{typeof(GeometricDistribution).Name}[p={P}]";
        }

    }
}
