using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DirectionalDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class SphereUniformDistribution :
        DirectionalDistribution<(ddouble x, ddouble y, ddouble z), (double x, double y, double z)> {

        private static readonly ddouble pdf = 1d / (4d * Pi);

        public SphereUniformDistribution() { }

        public override ddouble PDF((ddouble x, ddouble y, ddouble z) v) {
            ddouble r = Hypot(v.x, v.y, v.z);

            if (!(IsFinite(r) && r > 0d)) {
                return NaN;
            }

            return pdf;
        }

        public override (double x, double y, double z) Sample(Random random) {
            double r = random.NextUniform();
            double theta = 2 * random.NextUniformOpenInterval1();

            double w = (2 * r - 1);
            double t = double.Sqrt(1 - w * w);

            (double s, double c) = double.SinCosPi(theta);

            double x = t * c;
            double y = t * s;
            double z = w;

            return (x, y, z);
        }

        public override (ddouble x, ddouble y, ddouble z) Mean => (NaN, NaN, NaN);

        public override (ddouble x, ddouble y, ddouble z) Mode => (NaN, NaN, NaN);

        public override ddouble Variance => 1d;

        public override ddouble Skewness => 0d;

        public override ddouble Entropy => Log(4d * Pi);

        public override string ToString() {
            return $"{typeof(SphereUniformDistribution).Name}[]";
        }

        public override string Formula => "p(x) := 1 / (4 * pi)";
    }
}