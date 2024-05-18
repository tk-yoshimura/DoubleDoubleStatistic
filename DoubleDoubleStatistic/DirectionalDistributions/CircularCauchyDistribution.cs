using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DirectionalDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class CircularCauchyDistribution :
        DirectionalDistribution<ddouble, double>,
        IFittableDirectionalDistribution<CircularCauchyDistribution, ddouble, double> {

        public ddouble Gamma { get; }
        public ddouble Mu { get; }

        private readonly ddouble pdf_norm, cosh_gamma;

        public CircularCauchyDistribution(ddouble gamma, ddouble mu) {
            ParamAssert.ValidateShape(nameof(gamma), gamma >= 0d && gamma < 256d);
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));

            static ddouble round_mu(ddouble mu) {
                mu += PI;

                if (mu >= 0d) {
                    mu %= 2d * PI;
                }
                else {
                    mu = 2d * PI - ((-mu) % (2d * PI));
                }

                mu -= PI;

                return mu;
            }

            Gamma = gamma;
            Mu = gamma > 0d ? round_mu(mu) : 0d;

            pdf_norm = Sinh(gamma) / (2d * PI);
            cosh_gamma = Cosh(gamma);
        }

        public override ddouble PDF(ddouble x) {
            if (!IsFinite(x)) {
                return NaN;
            }

            if (IsFinite(cosh_gamma)) {
                ddouble pdf = pdf_norm / (cosh_gamma - Cos(x - Mu));

                return pdf;
            }
            else {
                return RcpPI * 0.5d;
            }
        }

        public override double Sample(Random random) {
            const double pi2 = 2 * double.Pi;

            double u = random.NextUniformOpenInterval01() - 0.5d;
            double r = double.TanPi(u);

            double v = r * (double)Gamma + (double)Mu + double.Pi;

            double theta = (v >= 0) ? (v % pi2) : (pi2 - ((-v) % pi2));

            theta -= double.Pi;

            if (double.IsFinite(theta)) {
                return theta;
            }
            else {
                return random.NextUniform() * pi2 - double.Pi;
            }
        }

        public override ddouble Mean => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance => 1d - Exp(-Gamma);

        public override ddouble Skewness => 0d;

        public override ddouble Entropy =>
            Log(2d * PI * (1d - Exp(-2d * Gamma)));

        public static CircularCauchyDistribution? Fit(IEnumerable<double> samples)
            => Fit(samples.Select(v => (ddouble)v));

        public static CircularCauchyDistribution? Fit(IEnumerable<ddouble> samples) {
            if (samples.Count() < 2) {
                return null;
            }

            ddouble x = samples.Select(Cos).Mean(), y = samples.Select(Sin).Mean();
            ddouble r = x * x + y * y;

            if (!(r > 0d)) {
                return null;
            }

            int n = samples.Count();

            ddouble re2 = (n * r - 1d) / (ddouble)(n - 1);

            ddouble gamma = Log(1d / re2) * 0.5d;
            ddouble mu = Atan2(y, x);

            try {
                return new CircularCauchyDistribution(gamma, mu);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(CircularCauchyDistribution).Name}[gamma={Gamma},mu={Mu}]";
        }

        public override string Formula => "p(x; gamma, mu) := sinh(gamma) / (cosh(gamma) - cos(x - mu)) / (2 * pi)";
    }
}