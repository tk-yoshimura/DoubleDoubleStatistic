using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
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
                ? kappa / (4d * Pi * Sinh(kappa))
                : (6d - kappa * kappa) / (24d * Pi);

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

            double w = (Kappa > 0d) ? (1 + double.Log(r + (1 - r) * exp_m2kappa) * kappa_inv) : (2 * r - 1);
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

        public override ddouble Variance => Kappa > 1e-15d
            ? 1d - BesselI(1.5, Kappa) / BesselI(0.5, Kappa)
            : 1d - Kappa / 3d + Cube(Kappa) / 45d;

        public override ddouble Skewness => 0d;

        public override ddouble Entropy => Kappa > 1e-15d
            ? -Log(pdf_norm) - Kappa * BesselI(1.5d, Kappa) / BesselI(0.5d, Kappa)
            : Log(4d * Pi);

        public static VonMisesFisherDistribution? Fit(IEnumerable<(double x, double y, double z)> samples)
            => Fit(samples.Select(v => ((ddouble)v.x, (ddouble)v.y, (ddouble)v.z)));

        public static VonMisesFisherDistribution? Fit(IEnumerable<(ddouble x, ddouble y, ddouble z)> samples) {
            if (samples.Count() < 2) {
                return null;
            }

            ddouble x = samples.Select(v => v.x).Mean(), y = samples.Select(v => v.y).Mean(), z = samples.Select(v => v.z).Mean();
            ddouble r = x * x + y * y + z * z;

            if (r < 1e-15d) {
                return new VonMisesFisherDistribution(kappa: 0d, (1d, 0d, 0d));
            }

            if (!(r > 0d)) {
                return null;
            }

            int n = samples.Count();

            ddouble re = Sqrt((n * r - 1d) / (ddouble)(n - 1));
            ddouble kappa = InverseKappaBp5.Value((double)re);

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


        internal static class InverseKappaBp5 {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0_0p25 = new(
                new ReadOnlyCollection<double>([
                    3.00000000000000000108e0,
                    -1.44398510234272622034e0,
                    -3.80785678971935392006e0,
                    -8.83103006427655465182e-2,
                    2.34571602840558150673e0,
                    7.57924083726435836333e-1,
                    -6.53592754377588964726e-1,
                    -4.47796684479304676600e-2,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.51867163255242473009e0,
                    -3.31327069891586297704e0,
                    -3.54849653613038491310e0,
                    1.24945211433061756267e0,
                    2.71503557272127255617e0,
                    2.08598203006309843510e-1,
                    -7.54638445859434347956e-1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p25_0p5 = new(
                new ReadOnlyCollection<double>([
                    6.42442885449375960513e-1,
                    1.96313307060568144083e0,
                    -1.24691438993448365947e0,
                    -9.17031469743059400653e0,
                    -5.19297470104047734130e0,
                    4.40223087372771053972e0,
                    1.52621810562575256173e0,
                    -3.64472106468736749202e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.08447469220094259324e0,
                    2.80155451319707694070e0,
                    -7.79259858788244224668e0,
                    -8.23697302719837043037e0,
                    3.00669531104333267150e0,
                    2.12649888665677904699e0,
                    -2.88958850643092893071e-1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p5_0p625 = new(
                new ReadOnlyCollection<double>([
                    7.21313727638361160765e-1,
                    3.07928348182007095301e0,
                    1.87547724374183695561e0,
                    -1.05755251261007640287e1,
                    -9.48937309185465452513e0,
                    6.30721387399584029053e0,
                    8.35223664587402199382e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.12066266977083830279e0,
                    6.99674402602881219365e0,
                    -7.28758175077822322330e0,
                    -1.33353788794498011496e1,
                    4.60120499566223307949e0,
                    1.50456784819099860689e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p625_0p75 = new(
                new ReadOnlyCollection<double>([
                    7.99561194094110575053e-1,
                    7.19705420416318655799e0,
                    2.80576842672937624678e1,
                    3.69515220876221229894e1,
                    -4.01414570473994211121e1,
                    -1.03332645251256247638e2,
                    -1.56413746602313712484e0,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.81866838365308647421e0,
                    4.26717748927997977071e1,
                    7.81314312060610423847e1,
                    8.45878417555142707152e0,
                    -1.02435418075931005427e2,
                    -1.50968955787410442842e1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p75_0p875 = new(
                new ReadOnlyCollection<double>([
                    8.88888711050418464523e-1,
                    1.13713927308750961561e1,
                    9.84311135959660034369e1,
                    4.91088335698386132666e2,
                    1.74305754459600536981e3,
                    2.56484886018081276536e3,
                    -8.74191379156823187162e2,
                    -7.20342562745210890423e3,
                    1.21490605001982350497e3,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.36817337021102964457e1,
                    1.22108179500487077607e2,
                    6.50972139033290361367e2,
                    2.45364311699828645450e3,
                    4.65469465453852923628e3,
                    1.86159217606757857711e3,
                    -7.14058800249783715166e3,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p875_0p9375 = new(
                new ReadOnlyCollection<double>([
                    8.88888711050418594571e-1,
                    8.35415041086725581985e0,
                    2.12984118098435786790e2,
                    2.59946261502963220403e3,
                    3.20277274650711873826e4,
                    3.02596156561378519798e5,
                    2.71989756607511669942e6,
                    5.90934951933613952309e2,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.50950677211332024271e0,
                    2.31254223902736093990e2,
                    2.71137389817281806660e3,
                    3.34324419346229394581e4,
                    3.08384226152427662560e5,
                    2.75735604640176521123e6,
                    -2.71953329740038908605e6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p9375_0p953125 = new(
                new ReadOnlyCollection<double>([
                    9.41176470588212887099e-1,
                    2.01584494484641404305e2,
                    2.61171387838187357026e4,
                    7.58058634103298750013e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.13242348919336338304e2,
                    2.75478754633235427197e4,
                    -2.61171387830527019713e4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p953125_0p96875 = new(
                new ReadOnlyCollection<double>([
                    9.69696969696969697020e-1,
                    6.42382878237546524999e-13,
                    -6.12039682646091082041e-13,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.69696969697632208022e-1,
                ])
            );

            public static double Value(double x) {
                if (!(x >= 0d && x <= 1d)) {
                    return double.NaN;
                }

                double v;
                if (x <= 0.25) {
                    v = x * ApproxUtil.Pade(x, pade_0_0p25);
                }
                else if (x <= 0.5) {
                    v = ApproxUtil.Pade(0.5 - x, pade_0p25_0p5);
                }
                else if (x <= 0.625) {
                    v = ApproxUtil.Pade(0.625 - x, pade_0p5_0p625);
                }
                else if (x <= 0.75) {
                    v = ApproxUtil.Pade(0.75 - x, pade_0p625_0p75);
                }
                else if (x <= 0.875) {
                    v = ApproxUtil.Pade(0.875 - x, pade_0p75_0p875);
                }
                else if (x <= 0.9375) {
                    v = ApproxUtil.Pade(x - 0.875, pade_0p875_0p9375);
                }
                else if (x <= 0.953125) {
                    v = ApproxUtil.Pade(x - 0.9375, pade_0p9375_0p953125);
                }
                else if (x <= 0.96875) {
                    v = ApproxUtil.Pade(0.96875 - x, pade_0p953125_0p96875);
                }
                else {
                    return 1d / (1d - x);
                }

                double y = v / (1d - v);

                return y;
            }
        }
    }
}