using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
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
            if (!IsFinite(x)) {
                return NaN;
            }

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
            ddouble kappa = InverseKappaB0.Value((double)re);
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

        internal static class InverseKappaB0 {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0_0p25 = new(
                new ReadOnlyCollection<double>([
                    2.00000000000000018653e0,
                    -1.72058848927204469017e-1,
                    -2.20193323789716439794e0,
                    1.52802474249326208734e-1,
                    4.42774151427558629882e-1,
                    -1.50795345807277563365e-2,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.91397057553642032899e0,
                    -1.77302546787772332123e0,
                    -2.08251728845593721498e0,
                    7.58006190610757013419e-1,
                    4.11372037325401515713e-1,
                    -4.64456526514403538606e-2,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p25_0p5 = new(
                new ReadOnlyCollection<double>([
                    5.36891226542936563377e-1,
                    1.46267169554112105143e0,
                    -1.73011417932378541414e0,
                    -6.65975232980687730742e0,
                    -1.78426237953357807897e0,
                    3.25542345451885021397e0,
                    5.37397300796598466139e-1,
                    -2.35763690126257963976e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.97770841055348325571e0,
                    2.21709319259326920200e0,
                    -6.43524638957145321520e0,
                    -4.85938428889713344436e0,
                    2.57753431679811561113e0,
                    1.09435897144106484291e0,
                    -2.06219028885077664720e-1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p5_0p625 = new(
                new ReadOnlyCollection<double>([
                    6.18691658052096569000e-1,
                    3.25473542236114638053e0,
                    2.87358138538520026961e0,
                    -8.70674503525510535109e0,
                    -1.21373064008980994018e1,
                    9.88150479456054440503e-1,
                    1.81378557910298019161e0,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.30897049687127087809e0,
                    1.11537367246605378725e1,
                    -1.41314847504653691723e0,
                    -1.40916402382485237299e1,
                    -2.20812551744684192668e0,
                    1.93438544839738598962e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p625_0p75 = new(
                new ReadOnlyCollection<double>([
                    7.03202549319484999654e-1,
                    7.23217783833818463351e0,
                    2.40136148210417755698e1,
                    1.05566686757383463281e1,
                    -5.65577774598028782868e1,
                    -3.38937995252136272244e1,
                    1.04277652289226937301e1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.13203215649701958361e1,
                    4.49148730586828132228e1,
                    5.45080231747120514405e1,
                    -3.35000351778295982191e1,
                    -4.91539385333586696877e1,
                    6.42547032835665244821e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p75_0p875 = new(
                new ReadOnlyCollection<double>([
                    8.12259862490935447322e-1,
                    3.09223255567435523568e1,
                    4.92062434651308923723e2,
                    4.61351848693759632089e3,
                    2.47064767823609056683e4,
                    7.01563689409377911901e4,
                    3.85028030771697592183e4,
                    -1.47350833692158224233e5,
                    -3.77418081391678006358e4,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.94198930330131393375e1,
                    6.55854926966429860722e2,
                    6.44851158993689396847e3,
                    3.73754132639570301637e4,
                    1.22028761348152527629e5,
                    1.46935873821197966461e5,
                    -1.01049810100707497899e5,
                    -6.50565284524637314063e4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p875_0p9375 = new(
                new ReadOnlyCollection<double>([
                    8.92222984417299119373e-1,
                    8.96579949791221509664e0,
                    2.68894986309391422815e2,
                    6.90134997042097273298e2,
                    1.99874038064822792979e4,
                    -2.67684976508270276839e4,
                    1.64557700535054543529e5,
                    -8.66470009781678974544e5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.17078478939650887753e1,
                    3.16884874897801270350e2,
                    1.25919756090508846306e3,
                    2.32918300879008741110e4,
                    5.31360437357427764368e3,
                    1.05042059307131505386e5,
                    -6.33803365061739785104e5,
                    -1.33822137537646350548e6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p9375_0p96875 = new(
                new ReadOnlyCollection<double>([
                    8.92222984417299170654e-1,
                    3.69214864038857978940e1,
                    3.98628942701323895272e2,
                    2.93345080569209757550e4,
                    2.43449158044635286402e5,
                    8.28975799005488194005e6,
                    2.57037977762187748091e7,
                    -7.59583976928501683020e6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.97224410620699470247e1,
                    3.76966550244120941349e2,
                    3.20912746947209625929e4,
                    2.17883543116474910518e5,
                    8.80070102012265449202e6,
                    1.31589538317983643115e7,
                    -6.68393224972722767274e7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_0p96875_1 = new(
                new ReadOnlyCollection<double>([
                    9.42071837304564167758e-1,
                    4.47103367913166423596e1,
                    9.71261841596315088301e1,
                    -5.77393480319984208322e3,
                    -2.11483419332225076062e4,
                    5.93350002185862925690e3,
                    -4.13221320110774087541e3,
                    5.20042325543523108460e3,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.56374232244610073566e1,
                    1.55714626367381271612e1,
                    -6.36520187796247990930e3,
                    -1.13287342088521732966e4,
                    5.35574193671963809645e4,
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
                    v = ApproxUtil.Pade(0.9375 - x, pade_0p875_0p9375);
                }
                else if (x <= 0.96875) {
                    v = ApproxUtil.Pade(x - 0.9375, pade_0p9375_0p96875);
                }
                else {
                    v = ApproxUtil.Pade(x - 0.96875, pade_0p96875_1);
                }

                double y = v / (1d - v);

                return y;
            }
        }
    }
}