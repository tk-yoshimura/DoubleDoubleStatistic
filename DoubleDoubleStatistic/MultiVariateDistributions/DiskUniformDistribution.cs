using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.MultiVariateDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class DiskUniformDistribution :
        MultiVariateDistribution {

        private static readonly ddouble pdf = PI;

        public DiskUniformDistribution() { }

        public override ddouble PDF(Vector x) {
            if (x.Dim != 2) {
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
            double theta = 2 * random.NextUniformOpenInterval1();
            double r = double.Sqrt(random.NextUniform());

            (double c, double s) = double.SinCosPi(2 * theta);

            return [r * c, r * s];
        }

        public override Vector Mean => Vector.Zero(2);

        public override Vector Mode => Vector.Invalid(2);

        public override Matrix Covariance => Matrix.Identity(2);

        public override ddouble Entropy => Log(pdf);

        public override string ToString() {
            return $"{typeof(DiskUniformDistribution).Name}[]";
        }

        public override string Formula => "p(x) := pi";
    }
}
