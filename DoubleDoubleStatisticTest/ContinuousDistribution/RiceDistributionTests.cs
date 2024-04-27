using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class RiceDistributionTests {
        readonly RiceDistribution dist_nu0 = new(nu: 0);
        readonly RiceDistribution dist_nu0p5 = new(nu: 0.5);
        readonly RiceDistribution dist_nu1 = new(nu: 1);
        readonly RiceDistribution dist_nu2 = new(nu: 2);
        readonly RiceDistribution dist_nu3 = new(nu: 3);
        readonly RiceDistribution dist_nu4 = new(nu: 4);
        readonly RiceDistribution dist_nu5 = new(nu: 5);
        readonly RiceDistribution dist_nu8 = new(nu: 8);
        readonly RiceDistribution dist_nu16 = new(nu: 16);
        readonly RiceDistribution dist_nu32 = new(nu: 32);
        readonly RiceDistribution dist_nu64 = new(nu: 64);
        readonly RiceDistribution dist_nu128 = new(nu: 128);


        RiceDistribution[] Dists => [
            dist_nu0,
            dist_nu0p5,
            dist_nu1,
            dist_nu2,
            dist_nu3,
            dist_nu4,
            dist_nu5,
            dist_nu8,
            dist_nu16,
            dist_nu32,
            dist_nu64,
            dist_nu128,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.Nu}");
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
        public void MeanTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (RiceDistribution dist in Dists) {
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
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);
                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);

                Console.WriteLine($"cdf(median)={cdf_median}");

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);
                ddouble ccdf_median = dist.CDF(dist.Median, Interval.Upper);

                Console.WriteLine($"ccdf(median)={ccdf_median}");

                Assert.IsTrue(ddouble.Abs(cdf_median + ccdf_median - 1) < 1e-28);

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
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
            foreach (RiceDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-280"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"quantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (RiceDistribution dist in Dists) {
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
            ddouble[] expected_nu1 = [
                0.000000000000000000e+00,
                7.552046500655910932e-02,
                1.492727699590049861e-01,
                2.195248696834362356e-01,
                2.846208141145958326e-01,
                3.430272019089722346e-01,
                3.933864868460277830e-01,
                4.345748244284223238e-01,
                4.657596075936404345e-01,
                4.864500452770387784e-01,
                4.965335208431001512e-01,
                4.962912184850060138e-01,
                4.863885328436609035e-01,
                4.678387335046347850e-01,
                4.419417186370438455e-01,
                4.102028731814614804e-01,
                3.742395128106319468e-01,
                3.356837647742091368e-01,
                2.960908444483334456e-01,
                2.568605934805660573e-01,
                2.191781089268603766e-01,
                1.839766990010553160e-01,
                1.519236825506342925e-01,
                1.234271041534655416e-01,
                9.865956527126350240e-02,
                7.759424150890376071e-02,
                6.004779506558975688e-02,
                4.572520923980331886e-02,
                3.426239575504532592e-02,
                2.526354513419463294e-02,
                1.833139806563943361e-02,
                1.308973917775960866e-02,
                9.198334505147100215e-03,
                6.361188960780916023e-03,
                4.329381375367986895e-03,
                2.899872640556384797e-03,
                1.911631533653251804e-03,
                1.240246133767455633e-03,
                7.919466046269983944e-04,
                4.977070524033340698e-04,
                3.078554166264388056e-04,
                1.874216941786176443e-04,
                1.123046259353706068e-04,
                6.623454745392799161e-05,
                3.844897956367901329e-05,
                2.196859745903573995e-05,
                1.235494604944120309e-05,
                6.839178331018433080e-06,
                3.726446699985958748e-06,
                1.998559098628154850e-06,
                1.055050368542850948e-06,
                5.482345194572054762e-07,
                2.804137278417618359e-07,
                1.411804181651745268e-07,
                6.996718776495825765e-08,
                3.413197959827889948e-08,
                1.638995916156048358e-08,
                7.747217086106660740e-09,
                3.604682890778634061e-09,
                1.650987640723501108e-09,
                7.443511092483323967e-10,
                3.303473891438315740e-10,
                1.443193687400246199e-10,
                6.206413016254036636e-11,
            ];
            ddouble[] expected_nu2 = [
                0.000000000000000000e+00,
                1.704855795299480004e-02,
                3.487466509371172502e-02,
                5.419468685585994427e-02,
                7.560500290056607064e-02,
                9.952822924512771952e-02,
                1.261675233458567635e-01,
                1.554724172002455718e-01,
                1.871197564053159734e-01,
                2.205129963687585681e-01,
                2.548021856313785749e-01,
                2.889254295894749247e-01,
                3.216705898121559049e-01,
                3.517536768581568940e-01,
                3.779081733742046811e-01,
                3.989777497815944485e-01,
                4.140038424479734469e-01,
                4.222995860945639168e-01,
                4.235027056998570960e-01,
                4.176020775780874272e-01,
                4.049354965639797577e-01,
                3.861593445551974035e-01,
                3.621938974529751443e-01,
                3.341505076522505169e-01,
                3.032485276951251141e-01,
                2.707304154818040987e-01,
                2.377829715727303728e-01,
                2.054712553302960043e-01,
                1.746896834695560730e-01,
                1.461324737291753162e-01,
                1.202833040070043352e-01,
                9.742210783993929257e-02,
                7.764552329091552185e-02,
                6.089674685576269147e-02,
                4.700040296358042047e-02,
                3.569842329336717773e-02,
                2.668368453027203727e-02,
                1.962910466118391803e-02,
                1.421088365596236483e-02,
                1.012546457045722902e-02,
                7.100500500102449077e-03,
                4.900602313346057864e-03,
                3.328912491851004646e-03,
                2.225631410897501838e-03,
                1.464560139236280442e-03,
                9.485693942304551565e-04,
                6.047048556209288876e-04,
                3.794324693318488519e-04,
                2.343398767501194389e-04,
                1.424565141801783401e-04,
                8.524057359035334589e-05,
                5.020438695340541661e-05,
                2.910525178415662427e-05,
                1.660880552243011786e-05,
                9.329225356765602213e-06,
                5.158168935424019273e-06,
                2.807313033018482630e-06,
                1.503951407009992381e-06,
                7.930977313031795642e-07,
                4.116913074229911864e-07,
                2.103637062065484374e-07,
                1.058098254993156319e-07,
                5.238887828796105394e-08,
                2.553359001146679840e-08,
            ];
            ddouble[] expected_nu4 = [
                0.000000000000000000e+00,
                4.424782578612365691e-05,
                1.029126484747513206e-04,
                1.930899558236616410e-04,
                3.374296077358279549e-04,
                5.673821711762056528e-04,
                9.269374683207367563e-04,
                1.476894478643834500e-03,
                2.299583626286775054e-03,
                3.503807458175033779e-03,
                5.229578102677784161e-03,
                7.652022685254657348e-03,
                1.098362483361466947e-02,
                1.547380652695265479e-02,
                2.140477401583982392e-02,
                2.908260259658111690e-02,
                3.882276164545776093e-02,
                5.092971244545044240e-02,
                6.567084922263435876e-02,
                8.324586393102649529e-02,
                1.037535208018973004e-01,
                1.271587099718456682e-01,
                1.532633633320175659e-01,
                1.816852033368381580e-01,
                2.118482182813110903e-01,
                2.429881257177692877e-01,
                2.741748857754380819e-01,
                3.043526386615202362e-01,
                3.323954419458754472e-01,
                3.571751286782213208e-01,
                3.776357959751964755e-01,
                3.928681596125969122e-01,
                4.021765094450081013e-01,
                4.051314039103230402e-01,
                4.016025440652498668e-01,
                3.917683201275386495e-01,
                3.761010548123698283e-01,
                3.553296267737998115e-01,
                3.303835727817596157e-01,
                3.023246100994989582e-01,
                2.722725589835474636e-01,
                2.413327793382136766e-01,
                2.105315069510403647e-01,
                1.807640504192388475e-01,
                1.527589449430216684e-01,
                1.270591472765473517e-01,
                1.040194783306041815e-01,
                8.381800154827519223e-02,
                6.647800793869625136e-02,
                5.189681099023422950e-02,
                3.987759940138137921e-02,
                3.016105059562954743e-02,
                2.245413417592777294e-02,
                1.645438412859328778e-02,
                1.186875950100858661e-02,
                8.426944990096659621e-03,
                5.889504033401333181e-03,
                4.051664638612741641e-03,
                2.743701126448972252e-03,
                1.828903859052690524e-03,
                1.200045811957310476e-03,
                7.751034257271091328e-04,
                4.928086001627274195e-04,
                3.084292162565575130e-04,
            ];

            foreach ((RiceDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 8, i++) {
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
            ddouble[] expected_nu1 = [
                0.000000000000000000e+00,
                4.729271894630231177e-03,
                1.880639342879497633e-02,
                4.190151794139462121e-02,
                7.347260204335158418e-02,
                1.127777465180364452e-01,
                1.588937819785770367e-01,
                2.107418397928363385e-01,
                2.671201962031782284e-01,
                3.267439485730068705e-01,
                3.882901995751652513e-01,
                4.504465499012426655e-01,
                5.119600008646963962e-01,
                5.716829852590024075e-01,
                6.286132627741157775e-01,
                6.819248492595081368e-01,
                7.309879399640865083e-01,
                7.753768203776015833e-01,
                8.148658846667304889e-01,
                8.494149447316090562e-01,
                8.791458768512627753e-01,
                9.043132230089514501e-01,
                9.252715976249541452e-01,
                9.424426575934428119e-01,
                9.562840284213597331e-01,
                9.672620255502674125e-01,
                9.758293633988949223e-01,
                9.824083968745416229e-01,
                9.873798657474628371e-01,
                9.910766616367291393e-01,
                9.937818341853046356e-01,
                9.957298970724771658e-01,
                9.971104672293480409e-01,
                9.980733421177564457e-01,
                9.987342560539930192e-01,
                9.991807246587152136e-01,
                9.994775594145243192e-01,
                9.996717923240204806e-01,
                9.997968819150645503e-01,
                9.998761713864379175e-01,
                9.999256378930536648e-01,
                9.999560126913207281e-01,
                9.999743707518463642e-01,
                9.999852915895663941e-01,
                9.999916860683603081e-01,
                9.999953714158154039e-01,
                9.999974620577125606e-01,
                9.999986294381463248e-01,
                9.999992710614921121e-01,
                9.999996181866988021e-01,
                9.999998030426687778e-01,
                9.999998999429587920e-01,
                9.999999499422363991e-01,
                9.999999753374615441e-01,
                9.999999880342130831e-01,
                9.999999942828857735e-01,
                9.999999973100767159e-01,
                9.999999987536899448e-01,
                9.999999994313675256e-01,
                9.999999997445216948e-01,
                9.999999998869706364e-01,
                9.999999999507562798e-01,
                9.999999999788727889e-01,
                9.999999999910728077e-01,
            ];
            ddouble[] expected_nu2 = [
                0.000000000000000000e+00,
                1.061426256896472707e-03,
                4.294622735221683557e-03,
                9.842506439092428214e-03,
                1.793063270833497938e-02,
                2.884876846538510578e-02,
                4.292621427237146292e-02,
                6.050209774785265254e-02,
                8.189230363059368800e-02,
                1.073551382725137504e-01,
                1.370581817146581394e-01,
                1.710489817828191039e-01,
                2.092322206032292409e-01,
                2.513556858405422423e-01,
                2.970067946574602713e-01,
                3.456205805530473407e-01,
                3.964990393880056252e-01,
                4.488406547089817678e-01,
                5.017779173567373796e-01,
                5.544198523112759958e-01,
                6.058960754724582731e-01,
                6.553987914358371247e-01,
                7.022194238874970296e-01,
                7.457772012322173572e-01,
                7.856379118373489900e-01,
                8.215220737635747783e-01,
                8.533027988740627956e-01,
                8.809945455712532292e-01,
                9.047346500041648643e-01,
                9.247599434794868678e-01,
                9.413808886923307195e-01,
                9.549555245145882054e-01,
                9.658651550685535270e-01,
                9.744932288392066377e-01,
                9.812083075940088994e-01,
                9.863514954448122651e-01,
                9.902282418510940376e-01,
                9.931040842498797261e-01,
                9.952036705168270103e-01,
                9.967122943602337859e-01,
                9.977791702628631132e-01,
                9.985217432121946679e-01,
                9.990304443966262493e-01,
                9.993734415260540072e-01,
                9.996010705621186654e-01,
                9.997497597011004711e-01,
                9.998453579426835880e-01,
                9.999058563575146108e-01,
                9.999435411055499490e-01,
                9.999666467780020662e-01,
                9.999805913954462611e-01,
                9.999888752602043018e-01,
                9.999937192253907092e-01,
                9.999965073596140952e-01,
                9.999980870568504487e-01,
                9.999989680738081699e-01,
                9.999994517419380502e-01,
                9.999997131180556886e-01,
                9.999998521587305511e-01,
                9.999999249664718315e-01,
                9.999999624961684042e-01,
                9.999999815392105473e-01,
                9.999999910510228851e-01,
                9.999999957279185336e-01,
            ];
            ddouble[] expected_nu4 = [
                0.000000000000000000e+00,
                2.692918193555151229e-06,
                1.165897776625299912e-05,
                2.972644863300309191e-05,
                6.217609133329218894e-05,
                1.176337304234636056e-04,
                2.093984807258677399e-04,
                3.572716383565128205e-04,
                5.899491443608236279e-04,
                9.480213445272867544e-04,
                1.487583084503094630e-03,
                2.284391629421370692e-03,
                3.438418398986394214e-03,
                5.078525658987331250e-03,
                7.366868690572534427e-03,
                1.050249099809482702e-02,
                1.472346410871511538e-02,
                2.030684843430904593e-02,
                2.756574354872311922e-02,
                3.684277919182856603e-02,
                4.849958951929212586e-02,
                6.290211753233294900e-02,
                8.040200248813046680e-02,
                1.013147793930555401e-01,
                1.258961166279701227e-01,
                1.543177709686642807e-01,
                1.866452776656954415e-01,
                2.228195498354046677e-01,
                2.626444853169678240e-01,
                3.057823467565866649e-01,
                3.517581034727018663e-01,
                3.999731507422663035e-01,
                4.497279363193708579e-01,
                5.002521313978687267e-01,
                5.507401980449212475e-01,
                6.003896318679546873e-01,
                6.484388752482379559e-01,
                6.942019413256422800e-01,
                7.370971546804989671e-01,
                7.766680495633153658e-01,
                8.125952832812128879e-01,
                8.446993126901355398e-01,
                8.729344337854908487e-01,
                8.973754995845264482e-01,
                9.181991381577009337e-01,
                9.356615539423276395e-01,
                9.500750109354305639e-01,
                9.617848976178757781e-01,
                9.711489152502291722e-01,
                9.785194801336716264e-01,
                9.842299537222252148e-01,
                9.885848703120820868e-01,
                9.918539632453258204e-01,
                9.942695217642349981e-01,
                9.960264487461801819e-01,
                9.972843267806140721e-01,
                9.981708185165082403e-01,
                9.987858037821054635e-01,
                9.992057667378002606e-01,
                9.994880699239202526e-01,
                9.996748718086483354e-01,
                9.997965490775954089e-01,
                9.998745684273372403e-01,
                9.999238133802261785e-01,
            ];

            foreach ((RiceDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 8, i++) {
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