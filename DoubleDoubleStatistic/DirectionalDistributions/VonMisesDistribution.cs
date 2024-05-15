using DoubleDouble;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DirectionalDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class VonMisesDistribution {

        public ddouble Mu { get; }
        public ddouble Kappa { get; }

        private readonly ddouble pdf_norm;

        private readonly double s;

        public VonMisesDistribution(ddouble mu, ddouble kappa) {
            if (!IsFinite(mu)) {
                throw new ArgumentOutOfRangeException(nameof(mu), "Invalid location parameter.");
            }
            if (!(kappa > 0d && kappa <= 256d)) {
                throw new ArgumentOutOfRangeException(nameof(kappa), "Invalid shape parameter.");
            }

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

            Mu = round_mu(mu);
            Kappa = kappa;

            pdf_norm = 1d / (2d * PI * BesselI(0, Kappa));

            s = (double)((kappa > 1.3d) ? (1d / Sqrt(kappa)) : (PI * Exp(-kappa)));
        }

        public ddouble PDF(ddouble x) {
            ddouble pdf = Exp(Kappa * Cos(x - Mu)) * pdf_norm;

            return pdf;
        }

        public double Sample(Random random) {
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

        public IEnumerable<double> Sample(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return Sample(random);
                count--;
            }
        }

        public virtual (ddouble min, ddouble max) Support => (-PI, PI);

        public virtual ddouble Mean => Mu;

        public virtual ddouble Median => Mu;

        public virtual ddouble Mode => Mu;

        public virtual ddouble Variance => 1d - BesselI(1, Kappa) / BesselI(0, Kappa);

        public virtual ddouble Skewness => 0d;

        public virtual ddouble Kurtosis => NaN;

        public virtual ddouble Entropy {
            get {
                ddouble b1 = BesselI(1, Kappa), b0 = BesselI(0, Kappa);

                return Log(2d * PI * b0) - Kappa * b1 / b0;
            }
        }

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
                return new VonMisesDistribution(mu, kappa);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(VonMisesDistribution).Name}[mu={Mu},kappa={Kappa}]";
        }

        public virtual string Formula => "p(x; mu, kappa) := exp(kappa * cos(x - mu)) / (2 * pi * bessel_i(0, kappa))";
    }
}