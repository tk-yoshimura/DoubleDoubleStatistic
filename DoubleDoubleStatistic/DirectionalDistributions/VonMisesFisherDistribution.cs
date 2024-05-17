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
    public class VonMisesFisherDistribution :
        DirectionalDistribution<(ddouble x, ddouble y, ddouble z), (double x, double y, double z)>,
        IFittableDirectionalDistribution<VonMisesFisherDistribution, (ddouble x, ddouble y, ddouble z), (double x, double y, double z)> {

        public (ddouble x, ddouble y, ddouble z) Mu { get; }
        public ddouble Kappa { get; }

        private readonly ddouble pdf_norm;

        private readonly double kappa_inv, exp_m2kappa;
        private readonly double qri, qrj, qij, qxx, qyy, qzz;

        public VonMisesFisherDistribution(ddouble kappa, (ddouble x, ddouble y, ddouble z) mu) {
            ddouble r = Hypot(mu.x, mu.y, mu.z);

            ParamAssert.ValidateShape(nameof(kappa), kappa >= 0d && kappa < 256d);
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(r) && r > 0d);

            Mu = kappa > 0d ? (mu.x / r, mu.y / r, mu.z / r) : (1d, 0d, 0d);
            Kappa = kappa;

            pdf_norm = kappa > 1e-15d
                ? kappa / (4d * PI * Sinh(kappa))
                : (6d - kappa * kappa) / (24d * PI);

            /* setup random gen. */
            {
                kappa_inv = (double)(1d / kappa);
                exp_m2kappa = (double)Exp(-2d * kappa);

                if (kappa > 0d) {
                    ddouble mx = Mu.x, my = Mu.y, mz = Mu.z;

                    ddouble angle = Acos(mz);
                    (ddouble s, ddouble c) = (Sin(angle * 0.5d), Cos(angle * 0.5d));
                    ddouble norm = Hypot(mx, my);

                    (ddouble qr, ddouble qi, ddouble qj) = (norm > 1e-30)
                        ? (c, s * -my / norm, s * mx / norm)
                        : (1d, 0d, 0d);

                    qri = (double)(qr * qi);
                    qij = (double)(qi * qj);
                    qrj = (double)(qr * qj);
                    qxx = (double)(qr * qr + qi * qi - qj * qj);
                    qyy = (double)(qr * qr - qi * qi + qj * qj);
                    qzz = (double)(qr * qr - qi * qi - qj * qj);
                }
                else {
                    (qri, qij, qrj, qxx, qyy, qzz) = (0d, 0d, 0.5d, 0d, 1d, 0d);
                }
            }
        }

        public override ddouble PDF((ddouble x, ddouble y, ddouble z) v) {
            ddouble r = Hypot(v.x, v.y, v.z);

            if (!(IsFinite(r) && r > 0d)) {
                return NaN;
            }

            ddouble pdf = Exp(Kappa * (v.x * Mu.x + v.y * Mu.y + v.z * Mu.z) / r) * pdf_norm;

            return pdf;
        }

        public override (double x, double y, double z) Sample(Random random) {
            double r = random.NextUniform();
            double theta = 2 * random.NextUniformOpenInterval1();

            double w = (Kappa > 0d) ? (1 + Math.Log(r + (1 - r) * exp_m2kappa) * kappa_inv) : (2 * r - 1);
            double t = double.Sqrt(1 - w * w);

            (double s, double c) = double.SinCosPi(theta);

            double x = t * c;
            double y = t * s;
            double z = w;

            (double, double, double) v = (
                x * qxx + 2 * (y * qij + z * qrj),
                y * qyy - 2 * (z * qri - x * qij),
                z * qzz - 2 * (x * qrj - y * qri)
            );

            return v;
        }

        public override (ddouble x, ddouble y, ddouble z) Mean => Mu;

        public override (ddouble x, ddouble y, ddouble z) Mode => Mu;

        public override ddouble Variance => 1d - BesselI(1.5, Kappa) / BesselI(0.5, Kappa);

        public override ddouble Skewness => 0d;

        public override ddouble Entropy => Kappa > 1e-15d
            ? -Log(pdf_norm) - Kappa * BesselI(1.5d, Kappa) / BesselI(0.5d, Kappa)
            : Log(4d * PI);

        public static VonMisesFisherDistribution? Fit(IEnumerable<(double x, double y, double z)> samples)
            => Fit(samples.Select(v => ((ddouble)v.x, (ddouble)v.y, (ddouble)v.z)));

        public static VonMisesFisherDistribution? Fit(IEnumerable<(ddouble x, ddouble y, ddouble z)> samples) {
            if (samples.Count() < 2) {
                return null;
            }

            ddouble x = samples.Select(v => v.x).Mean(), y = samples.Select(v => v.y).Mean(), z = samples.Select(v => v.z).Mean();
            ddouble r = x * x + y * y + z * z;

            if (!(r > 0d)) {
                return null;
            }

            int n = samples.Count();

            ddouble re = Sqrt((n * r - 1d) / (ddouble)(n - 1));

            ddouble t = BisectionMinimizeSearch1D.Search(
                t => {
                    ddouble kappa = t / (1d - t);

                    ddouble c = BesselI(1.5d, kappa) / BesselI(0.5d, kappa);
                    return Square(re - c);
                }, (0d, 257d / 256d), iter: 64
            );

            ddouble kappa = t / (1d - t);

            try {
                return new VonMisesFisherDistribution(kappa, (x, y, z));
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(VonMisesFisherDistribution).Name}[kappa={Kappa},mu={Mu}]";
        }

        public override string Formula => "p(x; kappa, mu) := exp(kappa * dot(mu, x)) / (kappa / (4 * pi * sinh(kappa)))";
    }
}