using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class NoncentralBetaDistributionTests {
        readonly NoncentralBetaDistribution dist_alpha_2_beta_2_lambda_1 = new(alpha: 2, beta: 2, lambda: 1);
        readonly NoncentralBetaDistribution dist_alpha_3_beta_2_lambda_2 = new(alpha: 3, beta: 2, lambda: 2);
        readonly NoncentralBetaDistribution dist_alpha_2_beta_3_lambda_3 = new(alpha: 2, beta: 3, lambda: 3);
        readonly NoncentralBetaDistribution dist_alpha_3_beta_3_lambda_4 = new(alpha: 3, beta: 3, lambda: 4);
        readonly NoncentralBetaDistribution dist_alpha_4_beta_5_lambda_6 = new(alpha: 4, beta: 5, lambda: 6);

        NoncentralBetaDistribution[] Dists => [
            dist_alpha_2_beta_2_lambda_1,
            dist_alpha_3_beta_2_lambda_2,
            dist_alpha_2_beta_3_lambda_3,
            dist_alpha_3_beta_3_lambda_4,
            dist_alpha_4_beta_5_lambda_6,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Beta={dist.Beta}");
                Console.WriteLine($"Lambda={dist.Lambda}");
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
            foreach (NoncentralBetaDistribution dist in Dists) {
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
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Console.WriteLine(dist.Mode);

                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode - 1e-15), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode + 1e-15), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
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
            foreach (NoncentralBetaDistribution dist in Dists) {
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
            foreach (NoncentralBetaDistribution dist in Dists) {
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
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
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

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-10"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-25);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
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

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-10"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"quantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-25);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (NoncentralBetaDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 40000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.1, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (NoncentralBetaDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha_2_beta_2_lambda_1 = [
                "0.8718640362787821151466971079602745644283",
                "1.472420230494374829026025036005887473795",
                "1.388381357484163275405606682365235093802",
            ];
            ddouble[] expected_dist_alpha_3_beta_2_lambda_2 = [
                "0.3113744366212743430841615448112361731005",
                "1.232015402541286641695217805450835296054",
                "2.032943840966665246893886839485401455297",
            ];
            ddouble[] expected_dist_alpha_2_beta_3_lambda_3 = [
                "0.8957417983286962858913583915654741379302",
                "1.717403003764450932299753231603307070301",
                "1.306569051768964436342637997852274912895",
            ];
            ddouble[] expected_dist_alpha_3_beta_3_lambda_4 = [
                "0.3681974493269622103906354252670802323613",
                "1.563487624978629866780976023186208686645",
                "1.995065326735427274466306429427776916314",
            ];
            ddouble[] expected_dist_alpha_4_beta_5_lambda_6 = [
                "0.3291288249475345417443404696716299337735",
                "2.104926093604821080007464367768625323786",
                "1.624711957250649682595390811994496344753",
            ];

            foreach ((NoncentralBetaDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_alpha_2_beta_2_lambda_1, expected_dist_alpha_2_beta_2_lambda_1),
                 (dist_alpha_3_beta_2_lambda_2, expected_dist_alpha_3_beta_2_lambda_2),
                 (dist_alpha_2_beta_3_lambda_3, expected_dist_alpha_2_beta_3_lambda_3),
                 (dist_alpha_3_beta_3_lambda_4, expected_dist_alpha_3_beta_3_lambda_4),
                 (dist_alpha_4_beta_5_lambda_6, expected_dist_alpha_4_beta_5_lambda_6),
            }) {
                for ((ddouble x, int i) = (1d / 4, 0); i < expecteds.Length; x += 1d / 4, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_alpha_2_beta_2_lambda_1 = [
                "0.1114160354290052587485386604475793515744",
                "0.4137379160066838362552467043322328438764",
                "0.7911446841529868943652060424651060311053",
            ];
            ddouble[] expected_dist_alpha_3_beta_2_lambda_2 = [
                "0.02537125039136309462167242216980442891930",
                "0.2084949142762177393638060901532182808707",
                "0.6365783744441072995324292123641156072143",
            ];
            ddouble[] expected_dist_alpha_2_beta_3_lambda_3 = [
                "0.1043171944829121186838784918044329102392",
                "0.4437662341180235822917976386791247385521",
                "0.8579580062544164362027435816280559819559",
            ];
            ddouble[] expected_dist_alpha_3_beta_3_lambda_4 = [
                "0.02857220922799107550671523606340931822363",
                "0.2586652320736703823718526508947771724228",
                "0.7536469220233026878494965169281281537580"
            ];
            ddouble[] expected_dist_alpha_4_beta_5_lambda_6 = [
                "0.01967540387750225203240559058246706169973",
                "0.3045364426366646496393769071648843991695",
                "0.8763868607925798439547571349019207118247",
            ];

            foreach ((NoncentralBetaDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_alpha_2_beta_2_lambda_1, expected_dist_alpha_2_beta_2_lambda_1),
                 (dist_alpha_3_beta_2_lambda_2, expected_dist_alpha_3_beta_2_lambda_2),
                 (dist_alpha_2_beta_3_lambda_3, expected_dist_alpha_2_beta_3_lambda_3),
                 (dist_alpha_3_beta_3_lambda_4, expected_dist_alpha_3_beta_3_lambda_4),
                 (dist_alpha_4_beta_5_lambda_6, expected_dist_alpha_4_beta_5_lambda_6),
            }) {
                for ((ddouble x, int i) = (1d / 4, 0); i < expecteds.Length; x += 1d / 4, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 2e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void PDFLargeLambdaTest() {
            NoncentralBetaDistribution dist_1_1_1024 = new(1, 1, 1024);
            NoncentralBetaDistribution dist_1_2_1024 = new(1, 2, 1024);
            NoncentralBetaDistribution dist_2_1_1024 = new(2, 1, 1024);
            NoncentralBetaDistribution dist_2_2_1024 = new(2, 2, 1024);
            NoncentralBetaDistribution dist_2_3_1024 = new(2, 3, 1024);

            ddouble expected_1_1_1024_0p25 = "2.195374722547311791531552508257025246077e-165";
            ddouble expected_1_1_1024_0p75 = "9.903006084673297083332221621893760960547e-54";
            ddouble expected_1_2_1024_0p25 = "2.156828026837469456587219435146930965594e-163";
            ddouble expected_1_2_1024_0p75 = "9.581094081687098867753253430730123510102e-52";
            ddouble expected_2_1_1024_0p25 = "5.530982828123072343005849342508009340892e-166";
            ddouble expected_2_1_1024_0p75 = "7.446546133799790923700462752047386488515e-54";
            ddouble expected_2_2_1024_0p25 = "5.475034809515519726613136328004947554097e-164";
            ddouble expected_2_2_1024_0p75 = "7.223053291934323105433442386807829565019e-52";
            ddouble expected_2_3_1024_0p25 = "2.750271638429302894165218767745937919358e-162";
            ddouble expected_2_3_1024_0p75 = "3.521098857077973767818983779892216024450e-50";

            Assert.IsTrue(ddouble.Abs(expected_1_1_1024_0p25 - dist_1_1_1024.PDF(0.25)) / expected_1_1_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_1_1024_0p75 - dist_1_1_1024.PDF(0.75)) / expected_1_1_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_2_1024_0p25 - dist_1_2_1024.PDF(0.25)) / expected_1_2_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_2_1024_0p75 - dist_1_2_1024.PDF(0.75)) / expected_1_2_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_1_1024_0p25 - dist_2_1_1024.PDF(0.25)) / expected_2_1_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_1_1024_0p75 - dist_2_1_1024.PDF(0.75)) / expected_2_1_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_2_1024_0p25 - dist_2_2_1024.PDF(0.25)) / expected_2_2_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_2_1024_0p75 - dist_2_2_1024.PDF(0.75)) / expected_2_2_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_3_1024_0p25 - dist_2_3_1024.PDF(0.25)) / expected_2_3_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_3_1024_0p75 - dist_2_3_1024.PDF(0.75)) / expected_2_3_1024_0p75 < 2e-29);
        }

        [TestMethod()]
        public void CDFLargeLambdaTest() {
            NoncentralBetaDistribution dist_1_1_1024 = new(1, 1, 1024);
            NoncentralBetaDistribution dist_1_2_1024 = new(1, 2, 1024);
            NoncentralBetaDistribution dist_2_1_1024 = new(2, 1, 1024);
            NoncentralBetaDistribution dist_2_2_1024 = new(2, 2, 1024);
            NoncentralBetaDistribution dist_2_3_1024 = new(2, 3, 1024);

            ddouble expected_1_1_1024_0p25 = "4.254602175479286417696807186544622569917e-168";
            ddouble expected_1_1_1024_0p75 = "1.9291570294818111201296535627065768104961e-56";
            ddouble expected_1_2_1024_0p25 = "4.158873626531002473298629024847368562094e-166";
            ddouble expected_1_2_1024_0p75 = "1.876105211171061314326088089732145948207481e-54";
            ddouble expected_2_1_1024_0p25 = "1.063650543869821604424201796636155642479e-168";
            ddouble expected_2_1_1024_0p75 = "1.4468677721113583400972401720299326078721e-56";
            ddouble expected_2_2_1024_0p25 = "1.047695785711774280357838769686613307842e-166";
            ddouble expected_2_2_1024_0p75 = "1.410696077808574381594809167729184292675291e-54";
            ddouble expected_2_3_1024_0p25 = "5.237614712491977171735599184473299662750e-165";
            ddouble expected_2_3_1024_0p75 = "6.9126820689692853495033331544052586464729012e-53";

            Assert.IsTrue(ddouble.Abs(expected_1_1_1024_0p25 - dist_1_1_1024.CDF(0.25)) / expected_1_1_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_1_1024_0p75 - dist_1_1_1024.CDF(0.75)) / expected_1_1_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_2_1024_0p25 - dist_1_2_1024.CDF(0.25)) / expected_1_2_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_2_1024_0p75 - dist_1_2_1024.CDF(0.75)) / expected_1_2_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_1_1024_0p25 - dist_2_1_1024.CDF(0.25)) / expected_2_1_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_1_1024_0p75 - dist_2_1_1024.CDF(0.75)) / expected_2_1_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_2_1024_0p25 - dist_2_2_1024.CDF(0.25)) / expected_2_2_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_2_1024_0p75 - dist_2_2_1024.CDF(0.75)) / expected_2_2_1024_0p75 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_3_1024_0p25 - dist_2_3_1024.CDF(0.25)) / expected_2_3_1024_0p25 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_3_1024_0p75 - dist_2_3_1024.CDF(0.75)) / expected_2_3_1024_0p75 < 2e-29);
        }
    }
}