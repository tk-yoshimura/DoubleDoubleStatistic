using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DirectionalDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class VonMisesDistribution :
        DirectionalDistribution<ddouble, double>,
        IFittableDirectionalDistribution<VonMisesDistribution, ddouble, double> {

        public ddouble Kappa { get; }
        public ddouble Mu { get; }

        private readonly ddouble pdf_norm;

        private readonly double s;

        public VonMisesDistribution(ddouble kappa, ddouble mu) {
            ParamAssert.ValidateShape(nameof(kappa), kappa >= 0d && kappa < 256d);
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

            Kappa = kappa;
            Mu = kappa > 0d ? round_mu(mu) : 0d;

            pdf_norm = 1d / (2d * PI * BesselI(0, kappa));

            s = (double)((kappa > 1.3d) ? (1d / Sqrt(kappa)) : (PI * Exp(-kappa)));
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = Exp(Kappa * Cos(x - Mu)) * pdf_norm;

            return pdf;
        }

        public override double Sample(Random random) {
            double kappa = (double)Kappa;

            double t;

            while (true) {
                double r1 = random.NextUniformOpenInterval01(), r2 = random.NextUniformOpenInterval01();
                t = s * (2 * r2 - 1) / r1;

                if (double.Abs(t) > double.Pi) {
                    continue;
                }

                if ((kappa * t * t < 4 - 4 * r1) || (kappa * double.Cos(t) >= 2 * double.Log(r1) + kappa)) {
                    break;
                }
            }

            double theta = t + (double)Mu;

            if (theta > double.Pi) {
                theta -= 2 * double.Pi;
            }
            else if (theta < -double.Pi) {
                theta += 2 * double.Pi;
            }

            return theta;
        }

        public override ddouble Mean => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance => 1d - BesselI(1, Kappa) / BesselI(0, Kappa);

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy =>
            -Log(pdf_norm) - Kappa * BesselI(1, Kappa) / BesselI(0, Kappa);

        public static VonMisesDistribution? Fit(IEnumerable<double> samples)
            => Fit(samples.Select(v => (ddouble)v));

        public static VonMisesDistribution? Fit(IEnumerable<ddouble> samples) {
            if (samples.Count() < 2) {
                return null;
            }

            ddouble x = samples.Select(Cos).Mean(), y = samples.Select(Sin).Mean();
            ddouble r = x * x + y * y;

            if (!(r > 0d)) {
                return null;
            }

            int n = samples.Count();

            ddouble re = Sqrt((n * r - 1d) / (ddouble)(n - 1));

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble kappa = t / (1d - t);

                    ddouble c = BesselI(1, kappa) / BesselI(0, kappa);
                    return Square(re - c);
                }, (0d, 257d / 256d), iter: 64
            );

            ddouble kappa = t / (1d - t);
            ddouble mu = Atan2(y, x);

            try {
                return new VonMisesDistribution(kappa, mu);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(VonMisesDistribution).Name}[kappa={Kappa},mu={Mu}]";
        }

        public override string Formula => "p(x; kappa, mu) := exp(kappa * cos(x - mu)) / (2 * pi * bessel_i(0, kappa))";
    }
}