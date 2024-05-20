using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.MultiVariateDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class BallUniformDistribution :
        MultiVariateDistribution {

        private static readonly ddouble pdf = 4d * PI / 3d;

        public BallUniformDistribution() { }

        public override ddouble PDF(Vector x) {
            if (x.Dim != 3) {
                throw new ArgumentException("Invalid argument dim.", nameof(x));
            }

            ddouble r = x.Norm;

            if (r <= 1d) {
                return pdf;
            }
            else if (r > 1d) {
                return 0d;
            }

            return NaN;
        }

        public override double[] Sample(Random random) {
            double theta = random.NextUniform() * 2 - 1;
            double t = double.Sqrt(1 - theta * theta);

            double phi = 2 * random.NextUniformOpenInterval1();
            double r = double.Cbrt(random.NextUniform());

            (double c, double s) = double.SinCosPi(2 * phi);

            return [r * t * c, r * t * s, r * theta];
        }

        public override Vector Mean => Vector.Zero(3);

        public override Vector Mode => Vector.Invalid(3);

        public override Matrix Covariance => Matrix.Identity(3) / 5d;

        public override ddouble Entropy => Log(pdf);

        public override string ToString() {
            return $"{typeof(BallUniformDistribution).Name}[]";
        }

        public override string Formula => "p(x) := 4 / 3 * pi";
    }
}
