using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.StableDistribution {
    [TestClass()]
    public class LevyDistributionTests {
        readonly LevyDistribution dist1 = new();
        readonly LevyDistribution dist2 = new(mu: 1, c: 3);

        LevyDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"C={dist.C}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode - 1e-4), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode + 1e-4), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (LevyDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.02, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (LevyDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 40000).ToArray();

                (LevyDistribution? dist_fit, ddouble error) = LevyDistribution.Fit(xs, (0.15, 0.8));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (LevyDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.NegativeInfinity)) && dist.PDF(ddouble.NegativeInfinity) >= 0d, "pdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MinValue)) && dist.PDF(ddouble.MinValue) >= 0d, "pdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MinValue / 2)) && dist.PDF(ddouble.MinValue / 2) >= 0d, "pdf(-lval / 2)");

                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.PositiveInfinity)) && dist.PDF(ddouble.PositiveInfinity) >= 0d, "pdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MaxValue)) && dist.PDF(ddouble.MaxValue) >= 0d, "pdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MaxValue / 2)) && dist.PDF(ddouble.MaxValue / 2) >= 0d, "pdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.PDF(ddouble.NaN)), "pdf(NaN)");

                Assert.IsTrue(dist.CDF(ddouble.NegativeInfinity, Interval.Lower) == 0d, "cdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue, Interval.Lower)) && dist.CDF(ddouble.MinValue, Interval.Lower) >= 0d, "cdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue / 2, Interval.Lower)) && dist.CDF(ddouble.MinValue / 2, Interval.Lower) >= 0d, "cdf(-lval / 2)");

                Assert.IsTrue(dist.CDF(ddouble.PositiveInfinity, Interval.Lower) == 1d, "cdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue, Interval.Lower)) && dist.CDF(ddouble.MaxValue, Interval.Lower) <= 1d, "cdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue / 2, Interval.Lower)) && dist.CDF(ddouble.MaxValue / 2, Interval.Lower) <= 1d, "cdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.CDF(ddouble.NaN, Interval.Lower)), "cdf(NaN)");

                Assert.IsTrue(dist.CDF(ddouble.NegativeInfinity, Interval.Upper) == 1d, "ccdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue, Interval.Upper)) && dist.CDF(ddouble.MinValue, Interval.Upper) <= 1d, "ccdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue / 2, Interval.Upper)) && dist.CDF(ddouble.MinValue / 2, Interval.Upper) <= 1d, "ccdf(-lval / 2)");

                Assert.IsTrue(dist.CDF(ddouble.PositiveInfinity, Interval.Upper) == 0d, "cdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue, Interval.Upper)) && dist.CDF(ddouble.MaxValue, Interval.Upper) >= 0d, "ccdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue / 2, Interval.Upper)) && dist.CDF(ddouble.MaxValue / 2, Interval.Upper) >= 0d, "ccdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.CDF(ddouble.NaN, Interval.Upper)), "ccdf(NaN)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(0d, Interval.Lower)), "quantile(0)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(1d, Interval.Lower)), "quantile(1)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(0d, Interval.Upper)), "cquantile(0)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(1d, Interval.Upper)), "cquantile(1)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.Epsilon, Interval.Lower)), "quantile(0+eps)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)), "quantile(1-eps)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.Epsilon, Interval.Upper)), "cquantile(0+eps)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)), "cquantile(1-eps)");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist1 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.653358828327364249e-01,
                4.319277321055045027e-01,
                4.579350180044978180e-01,
                4.151074974205947177e-01,
                3.627892871839331623e-01,
                3.153468637585579160e-01,
                2.752515104888593189e-01,
                2.419707245191433653e-01,
                2.143687681939464851e-01,
                1.913494655476511430e-01,
                1.720009130901971328e-01,
                1.555995547570865056e-01,
                1.415796865516508463e-01,
                1.294997276644349049e-01,
                1.190135187948509615e-01,
                1.098478223669305948e-01,
                1.017852217945228038e-01,
                9.465126089551177679e-02,
                8.830478300042056417e-02,
                8.263064759476451693e-02,
                7.753420634958216318e-02,
                7.293708326511037787e-02,
                6.877392602586679393e-02,
                6.498988524091373065e-02,
                6.153864312571569184e-02,
                5.838086006007908724e-02,
                5.548294138878026582e-02,
                5.281605149913378566e-02,
                5.035532024545615742e-02,
                4.807920006333168267e-02,
                4.596894195226054053e-02,
                4.400816584553744054e-02,
                4.218250640309396254e-02,
                4.047931943856413611e-02,
                3.888743737381713389e-02,
                3.739696455523393304e-02,
                3.599910515072166084e-02,
                3.468601781072824103e-02,
                3.345069242097598733e-02,
                3.228684517430722989e-02,
                3.118882890021845711e-02,
                3.015155615585357954e-02,
                2.917043303371039467e-02,
                2.824130200374549288e-02,
                2.736039239989285221e-02,
                2.652427739787252201e-02,
                2.572983652392880211e-02,
                2.497422289167438442e-02,
                2.425483449348873541e-02,
                2.356928897941269058e-02,
                2.291540144454115233e-02,
                2.229116481899731314e-02,
                2.169473251543559134e-02,
                2.112440303987983201e-02,
                2.057860631434220433e-02,
                2.005589149552707701e-02,
                1.955491610417351908e-02,
                1.907443630518219913e-02,
                1.861329820038600674e-02,
                1.817043001429727300e-02,
                1.774483506892328769e-02,
                1.733558545721956487e-02,
                1.694181633630587741e-02,
            ];
            ddouble[] expected_dist2 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                9.606648286403565978e-05,
                1.370231000044104010e-02,
                5.511196094424547498e-02,
                9.730434665928293370e-02,
                1.268656095908776427e-01,
                1.439759107018348250e-01,
                1.520385739690275073e-01,
                1.541803298037692782e-01,
                1.526450060014992727e-01,
                1.489197489412983633e-01,
                1.439596716727885306e-01,
                1.383691658068648966e-01,
                1.325267144191846769e-01,
                1.266663356887718639e-01,
                1.209297623946443828e-01,
                1.153997421040914434e-01,
                1.101213383329825851e-01,
                1.051156212528526340e-01,
                1.003885128144571676e-01,
                9.593652667283197877e-02,
                9.175050349628643964e-02,
                8.781804269015136555e-02,
                8.412508078846246684e-02,
                8.065690817304778382e-02,
                7.739881460743398833e-02,
                7.433648894461045153e-02,
                7.145625606464882373e-02,
                6.874520638433806452e-02,
                6.619125497909300604e-02,
                6.378315518255038563e-02,
                6.151048337760250256e-02,
                5.936360620416363670e-02,
                5.733363769673237992e-02,
                5.541239134535209249e-02,
                5.359233036000827372e-02,
                5.186651825236216623e-02,
                5.022857105612171802e-02,
                4.867261197056067107e-02,
                4.719322885055027977e-02,
                4.578543472609116621e-02,
                4.444463137700557787e-02,
                4.316657588814496599e-02,
                4.194735004917032190e-02,
                4.078333242818957954e-02,
                3.967117293161698949e-02,
                3.860776965754345730e-02,
                3.759024785250928563e-02,
                3.661594078897686261e-02,
                3.568237239103985020e-02,
                3.478724144763452319e-02,
                3.392840726484093461e-02,
                3.310387662118210311e-02,
                3.231179190177263516e-02,
                3.155042029850392560e-02,
                3.081814397405312558e-02,
            ];

            foreach ((LevyDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (-1, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist1 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                4.677734981047266226e-03,
                4.550026389635838941e-02,
                1.024704348597494097e-01,
                1.572992070502851614e-01,
                2.059032107320684091e-01,
                2.482130789899237300e-01,
                2.850494074026128644e-01,
                3.173105078629141484e-01,
                3.457785861511604164e-01,
                3.710933695226975582e-01,
                3.937686346429927209e-01,
                4.142161782425253236e-01,
                4.327675806677847570e-01,
                4.496917979688909162e-01,
                4.652088184521417924e-01,
                4.795001221869534813e-01,
                4.927166772270877848e-01,
                5.049850750938458255e-01,
                5.164122683960383764e-01,
                5.270892568655380916e-01,
                5.370939784426416175e-01,
                5.464935954065819335e-01,
                5.553463166714220911e-01,
                5.637028616507731016e-01,
                5.716076449533316062e-01,
                5.790997419539187785e-01,
                5.862136810731399805e-01,
                5.929800980174266822e-01,
                5.994262792958635622e-01,
                6.055766163353464293e-01,
                6.114529869535042517e-01,
                6.170750774519737636e-01,
                6.224606558934542289e-01,
                6.276258050283592960e-01,
                6.325851216960416412e-01,
                6.373518882339370695e-01,
                6.419382204050629870e-01,
                6.463551955394901682e-01,
                6.506129639327535852e-01,
                6.547208460185769408e-01,
                6.586874174078845012e-01,
                6.625205835400576060e-01,
                6.662276454096187628e-01,
                6.698153575994165720e-01,
                6.732899796599957076e-01,
                6.766573217164244536e-01,
                6.799227850521523120e-01,
                6.830913983096087438e-01,
                6.861678498552392647e-01,
                6.891565167793516355e-01,
                6.920614909359381617e-01,
                6.948866023724733498e-01,
                6.976354404528670727e-01,
                7.003113729368903861e-01,
                7.029175632453666944e-01,
                7.054569861112733875e-01,
                7.079324417918871903e-01,
                7.103465689955668072e-01,
                7.127018566581784231e-01,
                7.150006546880892655e-01,
                7.172451837847029221e-01,
                7.194375444233913619e-01,
                7.215797250891102799e-01,
            ];
            ddouble[] expected_dist2 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                9.633570086430964949e-07,
                5.320055051392504151e-04,
                4.677734981047266226e-03,
                1.430587843542964824e-02,
                2.845973691631056499e-02,
                4.550026389635838941e-02,
                6.407750645105955711e-02,
                8.326451666355044745e-02,
                1.024704348597494097e-01,
                1.213352503584822217e-01,
                1.396493986104293816e-01,
                1.572992070502851614e-01,
                1.742313882480249765e-01,
                1.904302638255243041e-01,
                2.059032107320684091e-01,
                2.206713619198469312e-01,
                2.347636628276692505e-01,
                2.482130789899237300e-01,
                2.610541895796574252e-01,
                2.733216782922981425e-01,
                2.850494074026128644e-01,
                2.962698714842864067e-01,
                3.070138980338251722e-01,
                3.173105078629141484e-01,
                3.271868777903057524e-01,
                3.366683676100388212e-01,
                3.457785861511604164e-01,
                3.545394797735015313e-01,
                3.629714323398304243e-01,
                3.710933695226975582e-01,
                3.789228628692691281e-01,
                3.864762307712327205e-01,
                3.937686346429927209e-01,
                4.008141693829344598e-01,
                4.076259477027809330e-01,
                4.142161782425253236e-01,
                4.205962375999267033e-01,
                4.267767365329836471e-01,
                4.327675806677847570e-01,
                4.385780260809999387e-01,
                4.442167301386065192e-01,
                4.496917979688909162e-01,
                4.550108249342736944e-01,
                4.601809354471203539e-01,
                4.652088184521417924e-01,
                4.701007598741285820e-01,
                4.748626723057679522e-01,
                4.795001221869534813e-01,
                4.840183547047819390e-01,
                4.884223166225936108e-01,
                4.927166772270877848e-01,
                4.969058475647679662e-01,
                5.009939981227101713e-01,
                5.049850750938458255e-01,
                5.088828153535194243e-01,

            ];

            foreach ((LevyDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (-1, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}