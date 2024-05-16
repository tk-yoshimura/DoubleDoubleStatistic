using DoubleDouble;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DiscreteDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class DiscreteUniformDistribution : DiscreteDistribution {

        public int A { get; }
        public int B { get; }
        public int N { get; }

        private readonly ddouble r, inv_r;

        public DiscreteUniformDistribution(int n) : this(0, n) { }

        public DiscreteUniformDistribution(int a, int b) {
            ParamAssert.ValidateShape(nameof(a), a < b);
            ParamAssert.ValidateShape(nameof(b), b - a >= 0);

            A = a;
            B = b;
            N = b - a;
            r = B - A;
            inv_r = 1d / r;
        }

        public override ddouble PMF(int k) {
            return k >= A && k < B ? inv_r : 0d;
        }

        public override int Sample(Random random) {
            return random.Next(A, B);
        }

        public override (int min, int max) Support => (A, B - 1);

        public override ddouble Mean => (A + B - 1) * 0.5d;

        public override ddouble Variance => checked((B - A - 1) * (B - A + 1)) / (ddouble)12;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => r > 1d
            ? (ddouble)checked(-6 * ((B - A) * (B - A) + 1)) / checked(5 * (B - A - 1) * (B - A + 1))
            : NegativeInfinity;

        public override ddouble Entropy => Log(r);

        public override string Formula => "f(k; a, b) := 1 / (b - a)";

        public override string ToString() {
            return $"{typeof(DiscreteUniformDistribution).Name}[a={A},b={B}]";
        }

    }
}
