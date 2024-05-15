using DoubleDouble;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.DirectionalDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class VonMisesFisherDistribution {

        public (ddouble x, ddouble y, ddouble z) Mu { get; }
        public ddouble Kappa { get; }

        private readonly ddouble pdf_norm;

        private readonly double kappa_inv, exp_m2kappa;
        private readonly double qri, qrj, qij, qxx, qyy, qzz;

        public VonMisesFisherDistribution((ddouble x, ddouble y, ddouble z) mu, ddouble kappa) {
            if (!(kappa > 0d && kappa <= 256d)) {
                throw new ArgumentOutOfRangeException(nameof(kappa), "Invalid shape parameter.");
            }

            ddouble r = Hypot(mu.x, mu.y, mu.z);

            if (!(IsFinite(r) && r > 0d)) {
                throw new ArgumentOutOfRangeException(nameof(mu), "Invalid location parameter.");
            }

            Mu = (mu.x / r, mu.y / r, mu.z / r);
            Kappa = kappa;

            pdf_norm = kappa / (4d * PI * Sinh(kappa));

            /* setup random gen. */
            {
                kappa_inv = (double)(1d / kappa);
                exp_m2kappa = (double)Exp(-2d * kappa);

                ddouble mx = Mu.x, my = Mu.y, mz = Mu.z;

                ddouble angle = Acos(mz);
                (ddouble s, ddouble c) = (Sin(angle * 0.5d), Cos(angle * 0.5d));
                ddouble norm = Hypot(mx, my);

                (ddouble qr, ddouble qi, ddouble qj) = (norm > 0d)
                    ? (c, s * -my / norm, s * mx / norm)
                    : (1d, 0d, 0d);

                qri = (double)(qr * qi);
                qij = (double)(qi * qj);
                qrj = (double)(qr * qj);
                qxx = (double)(qr * qr + qi * qi - qj * qj);
                qyy = (double)(qr * qr - qi * qi + qj * qj);
                qzz = (double)(qr * qr - qi * qi - qj * qj);
            }
        }

        public ddouble PDF((ddouble x, ddouble y, ddouble z) v) {
            ddouble r = Hypot(v.x, v.y, v.z);

            if (!(IsFinite(r) && r > 0d)) {
                return NaN;
            }

            ddouble pdf = Exp(Kappa * (v.x * Mu.x + v.y * Mu.y + v.z * Mu.z) / r) * pdf_norm;

            return pdf;
        }

        public (double x, double y, double z) Sample(Random random) {
            double r = random.NextUniform();
            double theta = 2 * random.NextUniformOpenInterval1();

            double w = 1 + double.Log(r + (1 - r) * exp_m2kappa) * kappa_inv;
            double t = double.Sqrt(1 - w * w);

            (double s, double c) = double.SinCosPi(theta);

            double x = t * c;
            double y = t * s;
            double z = w;

            (double, double, double) v = (x * qxx + 2 * (y * qij + z * qrj), y * qyy - 2 * (z * qri - x * qij), z * qzz - 2 * (x * qrj - y * qri));

            return v;
        }

        public IEnumerable<(double x, double y, double z)> Sample(Random random, int count) {
            ArgumentOutOfRangeException.ThrowIfNegative(count, nameof(count));

            while (count > 0) {
                yield return Sample(random);
                count--;
            }
        }

        public virtual (ddouble x, ddouble y, ddouble z) Mean => Mu;

        public virtual (ddouble x, ddouble y, ddouble z) Mode => Mu;

        public virtual ddouble Entropy => -Log(pdf_norm) - Kappa * BesselI(1.5d, Kappa) / BesselI(0.5d, Kappa);

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

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble kappa = t / (1d - t);

                    ddouble c = BesselI(1.5d, kappa) / BesselI(0.5d, kappa);
                    return Square(re - c);
                }, (0d, 257d / 256d), iter: 64
            );

            ddouble kappa = t / (1d - t);

            try {
                return new VonMisesFisherDistribution((x, y, z), kappa);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(VonMisesFisherDistribution).Name}[mu={Mu},kappa={Kappa}]";
        }

        public virtual string Formula => "p(x; mu, kappa) := exp(kappa * dot(mu, x)) / (kappa / (4 * pi * sinh(kappa)))";
    }
}