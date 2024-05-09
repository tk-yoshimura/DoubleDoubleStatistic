using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LogarithmicDistribution : DiscreteDistribution {
        const int random_gen_max_index = 65536;

        public ddouble P { get; }

        private readonly ddouble p, q, lnq;

        private readonly double fr;

        public LogarithmicDistribution(ddouble p) {
            ValidateShape(p, p => p > 0d && p < 1d);

            P = p;

            this.p = p;
            q = 1d - p;
            lnq = Log(q);

            fr = (double)(-p / lnq);
        }

        public override ddouble PMF(int k) {
            return k >= 1 ? -Pow(p, k) / (lnq * k) : 0d;
        }

        public override int Sample(Random random) {
            int k = 1;
            double p = (double)P, f = fr, r = random.NextUniformOpenInterval1();

            while(true) {
                r -= f;
                if (r <= 0 || k >= random_gen_max_index || f < double.Epsilon) {
                    return k;
                }

                f *= p * k / (k + 1);
                k++;
            }
        }

        public override (int min, int max) Support => (1, int.MaxValue);

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
