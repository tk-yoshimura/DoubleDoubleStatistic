using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class ZipfDistribution : DiscreteDistribution {
        const int random_gen_max_index = 65536;

        private Roulette roulette = null;

        public ddouble S { get; }

        private readonly ddouble pdfnorm;

        public ZipfDistribution(ddouble s) {
            ValidateShape(s, s => s > 0d);

            S = s;

            pdfnorm = 1d / RiemannZeta(s + 1d);
        }

        public override ddouble PMF(int k) {
            return k >= 1 ? pdfnorm / Pow(k, S + 1d) : 0d;
        }

        public override int Sample(Random random) {
            if (roulette is null) {
                SetupRoulette();
            }

            return roulette.NextIndex(random);
        }

        public void SetupRoulette() {
            List<double> probs = [];
            ddouble ps = 0d;
            for (int k = 0; k <= random_gen_max_index && (double)ps < 1d; k++) {
                ddouble p = PMF(k);
                ps += p;

                probs.Add((double)p);
            }

            roulette = new(probs);
        }

        public override (int min, int max) Support => (1, int.MaxValue);

        public override ddouble Mean => S > 1d
            ? RiemannZeta(S) / RiemannZeta(S + 1d)
            : PositiveInfinity;

        public override ddouble Variance => S > 2d
            ? RiemannZeta(S - 1d) / RiemannZeta(S + 1d) - Square(Mean)
            : PositiveInfinity;

        public override ddouble Skewness {
            get {
                if (S <= 3d) {
                    return PositiveInfinity;
                }

                ddouble sm2 = RiemannZeta(S - 2d), sm1 = RiemannZeta(S - 1d);
                ddouble s0 = RiemannZeta(S), sp1_inv = 1d / RiemannZeta(S + 1d);

                return sp1_inv * (sm2 + sp1_inv * s0 * (-3d * sm1 + sp1_inv * 2d * s0 * s0)) / ExMath.Pow3d2(Variance);
            }
        }

        public override ddouble Kurtosis {
            get {
                if (S <= 4d) {
                    return PositiveInfinity;
                }

                ddouble sm3 = RiemannZeta(S - 3d), sm2 = RiemannZeta(S - 2d), sm1 = RiemannZeta(S - 1d);
                ddouble s0 = RiemannZeta(S), sp1 = RiemannZeta(S + 1d);

                return (-3d * Pow(s0, 4) + 6d * sm1 * sp1 * s0 * s0 - 4d * sm2 * sp1 * sp1 * s0 + sm3 * Cube(sp1)) / Square(s0 * s0 - sm1 * sp1) - 3d;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, 1, 65536);

        public override string Formula => "f(k; s) := k^(- s - 1) / zeta(s + 1)";

        public override string ToString() {
            return $"{typeof(ZipfDistribution).Name}[s={S}]";
        }

    }
}
