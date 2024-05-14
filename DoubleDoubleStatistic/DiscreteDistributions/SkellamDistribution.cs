using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class SkellamDistribution : DiscreteDistribution {
        public ddouble Mu1 { get; }
        public ddouble Mu2 { get; }

        private readonly ddouble mu12, lnmud;

        private readonly int mean_n;

        private readonly ReadOnlyCollection<ddouble> besseli_table;

        private readonly PoissonDistribution randam_gen_poisson_dist_mu1, randam_gen_poisson_dist_mu2;

        public SkellamDistribution(ddouble mu1, ddouble mu2) {
            ValidateShape(mu1, mu => mu > 0d && mu <= 256d);
            ValidateShape(mu2, mu => mu > 0d && mu <= 256d);

            Mu1 = mu1;
            Mu2 = mu2;

            mu12 = mu1 + mu2;
            lnmud = Log(mu1 / mu2);

            mean_n = (int)Truncate(mu1 - mu2);

            ddouble x = 2d * Sqrt(mu1 * mu2);

            besseli_table = ExMath.BesselIIntegerNuTable(x);

            randam_gen_poisson_dist_mu1 = new(mu1);
            randam_gen_poisson_dist_mu2 = new(mu2);
        }

        public override ddouble PMF(int k) {
            if (k <= -besseli_table.Count || k >= besseli_table.Count) {
                return 0d;
            }

            ddouble bi = besseli_table[int.Abs(k)];

            return Exp(-mu12 + lnmud * k * 0.5d) * bi;
        }

        public override int Sample(Random random) {
            long n1 = randam_gen_poisson_dist_mu1.Sample(random);
            long n2 = randam_gen_poisson_dist_mu2.Sample(random);

            int n = (int)long.Clamp(n1 - n2, int.MinValue, int.MaxValue);

            return n;
        }

        public override (int min, int max) Support => (int.MinValue, int.MaxValue);

        public override ddouble Mean => Mu1 - Mu2;

        public override ddouble Variance => Mu1 + Mu2;

        public override ddouble Skewness => Mean / ExMath.Pow3d2(Variance);

        public override ddouble Kurtosis => 1d / (Mu1 + Mu2);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            DiscreteEntropy.Sum(this, mean_n, mean_n + 65535) + DiscreteEntropy.Sum(this, mean_n - 1, mean_n - 65536);

        public override string Formula => "f(k; mu1, mu2) := exp(- mu1 - mu2) * (mu1 / mu2)^(k / 2) * bessel_i(2 * sqrt(mu1 * mu2))";

        public override string ToString() {
            return $"{typeof(SkellamDistribution).Name}[mu1={Mu1},mu2={Mu1}]";
        }

    }
}
