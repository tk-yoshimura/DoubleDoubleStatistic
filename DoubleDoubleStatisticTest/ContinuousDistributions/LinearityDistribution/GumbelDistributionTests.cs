using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class GumbelDistributionTests {
        readonly GumbelDistribution dist_a1b1 = new(mu: 1, sigma: 1);
        readonly GumbelDistribution dist_a2b1 = new(mu: 2, sigma: 1);
        readonly GumbelDistribution dist_a1b2 = new(mu: 1, sigma: 2);
        readonly GumbelDistribution dist_a2b2 = new(mu: 2, sigma: 2);
        readonly GumbelDistribution dist_a3b4 = new(mu: 3, sigma: 4);

        GumbelDistribution[] Dists => [
            dist_a1b1,
            dist_a2b1,
            dist_a1b2,
            dist_a2b2,
            dist_a3b4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Sigma={dist.Sigma}");
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
            foreach (GumbelDistribution dist in Dists) {
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
            foreach (GumbelDistribution dist in Dists) {
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
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (GumbelDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (GumbelDistribution dist in Dists) {
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
            foreach (GumbelDistribution dist in Dists) {
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

            foreach (GumbelDistribution dist in Dists) {

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
        public void IrregularValueTest() {
            foreach (GumbelDistribution dist in Dists) {
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
            ddouble[] expected_dist_a1b1 = [
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
                9.027618899817316972e-287,
                5.549184995570460952e-223,
                2.505347870200592702e-173,
                1.111535258546083656e-134,
                1.319057421837700937e-104,
                3.297366349356508657e-81,
                5.205427108495642505e-63,
                7.333002349558938559e-49,
                7.250730578601260254e-38,
                2.508174882337907440e-29,
                1.060480399704279740e-22,
                1.451887713099222174e-17,
                1.374588275433549892e-13,
                1.625005195077556459e-10,
                3.800542504044356834e-08,
                2.516524225309105680e-06,
                6.236577187661992917e-05,
                7.189377113313632833e-04,
                4.566281420127915255e-03,
                1.823153164476988553e-02,
                5.070711360998073070e-02,
                1.064219999008989004e-01,
                1.793740787340172005e-01,
                2.548704208230366808e-01,
                3.170419210779421570e-01,
                3.555727473819441653e-01,
                3.678794411714423340e-01,
                3.574353461721825331e-01,
                3.307042988904180802e-01,
                2.945323152403546696e-01,
                2.546463800435824765e-01,
                2.151317179402431057e-01,
                1.785065185131209375e-01,
                1.460547184694531153e-01,
                1.182049515931431344e-01,
                9.485563027712225204e-02,
                7.561617991742651534e-02,
                5.996897935496824095e-02,
                4.736900967790790701e-02,
                3.729954287267522872e-02,
                2.929913213328151633e-02,
                2.297111444931936378e-02,
                1.798322969671364152e-02,
                1.406220982483213922e-02,
                1.098626968240280256e-02,
                8.577166239058568731e-03,
                6.692699677535513675e-03,
                5.220054072633552832e-03,
                4.070103819247028118e-03,
                3.172666806755189330e-03,
                2.472615573014907715e-03,
                1.926731077808604715e-03,
                1.501180561853306247e-03,
                1.169509464005654429e-03,
            ];
            ddouble[] expected_dist_a2b1 = [
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
                9.027618899817316972e-287,
                5.549184995570460952e-223,
                2.505347870200592702e-173,
                1.111535258546083656e-134,
                1.319057421837700937e-104,
                3.297366349356508657e-81,
                5.205427108495642505e-63,
                7.333002349558938559e-49,
                7.250730578601260254e-38,
                2.508174882337907440e-29,
                1.060480399704279740e-22,
                1.451887713099222174e-17,
                1.374588275433549892e-13,
                1.625005195077556459e-10,
                3.800542504044356834e-08,
                2.516524225309105680e-06,
                6.236577187661992917e-05,
                7.189377113313632833e-04,
                4.566281420127915255e-03,
                1.823153164476988553e-02,
                5.070711360998073070e-02,
                1.064219999008989004e-01,
                1.793740787340172005e-01,
                2.548704208230366808e-01,
                3.170419210779421570e-01,
                3.555727473819441653e-01,
                3.678794411714423340e-01,
                3.574353461721825331e-01,
                3.307042988904180802e-01,
                2.945323152403546696e-01,
                2.546463800435824765e-01,
                2.151317179402431057e-01,
                1.785065185131209375e-01,
                1.460547184694531153e-01,
                1.182049515931431344e-01,
                9.485563027712225204e-02,
                7.561617991742651534e-02,
                5.996897935496824095e-02,
                4.736900967790790701e-02,
                3.729954287267522872e-02,
                2.929913213328151633e-02,
                2.297111444931936378e-02,
                1.798322969671364152e-02,
                1.406220982483213922e-02,
                1.098626968240280256e-02,
                8.577166239058568731e-03,
                6.692699677535513675e-03,
                5.220054072633552832e-03,
                4.070103819247028118e-03,
                3.172666806755189330e-03,
            ];
            ddouble[] expected_dist_a1b2 = [
                3.625365289300630127e-38,
                1.255233154720231378e-33,
                1.254087441168953720e-29,
                4.183950880775744723e-26,
                5.302401998521398699e-23,
                2.860099425024682965e-20,
                7.259438565496110871e-18,
                9.473431039362024917e-16,
                6.872941377167749459e-14,
                2.970083759098640948e-12,
                8.125025975387782294e-11,
                1.484716081193370294e-09,
                1.900271252022178417e-08,
                1.776285433889975965e-07,
                1.258262112654552840e-06,
                6.978187351855016993e-06,
                3.118288593830996459e-05,
                1.151631929225159430e-04,
                3.594688556656816416e-04,
                9.672586044006541194e-04,
                2.283140710063957628e-03,
                4.800830760739250962e-03,
                9.115765822384942763e-03,
                1.581864110300537640e-02,
                2.535355680499036535e-02,
                3.788409532148280523e-02,
                5.321099995044945019e-02,
                7.076689095935798091e-02,
                8.968703936700860024e-02,
                1.089329771368782868e-01,
                1.274352104115183404e-01,
                1.442232059743190353e-01,
                1.585209605389710785e-01,
                1.697992206073497823e-01,
                1.777863736909720827e-01,
                1.824469864020708543e-01,
                1.839397205857211670e-01,
                1.825658985584724392e-01,
                1.787176730860912666e-01,
                1.728317917757161704e-01,
                1.653521494452090401e-01,
                1.567021186831329960e-01,
                1.472661576201773348e-01,
                1.373793501883596724e-01,
                1.273231900217912382e-01,
                1.173259115464438673e-01,
                1.075658589701215528e-01,
                9.817665938293082539e-02,
                8.925325925656046877e-02,
                8.085815309826707709e-02,
                7.302735923472655766e-02,
                6.577587511923496666e-02,
                5.910247579657156719e-02,
                5.299381084138009290e-02,
                4.742781513856112602e-02,
                4.237648499479453046e-02,
                3.780808995871325767e-02,
                3.368889781496139935e-02,
                2.998448967748412047e-02,
                2.666073700370180655e-02,
                2.368450483895395350e-02,
                2.102413712941540030e-02,
                1.864977143633761436e-02,
                1.653352239387048325e-02,
            ];
            ddouble[] expected_dist_a2b2 = [
                2.602713554247821253e-63,
                8.606044595302582537e-56,
                3.666501174779469279e-49,
                2.560592649691453994e-43,
                3.625365289300630127e-38,
                1.255233154720231378e-33,
                1.254087441168953720e-29,
                4.183950880775744723e-26,
                5.302401998521398699e-23,
                2.860099425024682965e-20,
                7.259438565496110871e-18,
                9.473431039362024917e-16,
                6.872941377167749459e-14,
                2.970083759098640948e-12,
                8.125025975387782294e-11,
                1.484716081193370294e-09,
                1.900271252022178417e-08,
                1.776285433889975965e-07,
                1.258262112654552840e-06,
                6.978187351855016993e-06,
                3.118288593830996459e-05,
                1.151631929225159430e-04,
                3.594688556656816416e-04,
                9.672586044006541194e-04,
                2.283140710063957628e-03,
                4.800830760739250962e-03,
                9.115765822384942763e-03,
                1.581864110300537640e-02,
                2.535355680499036535e-02,
                3.788409532148280523e-02,
                5.321099995044945019e-02,
                7.076689095935798091e-02,
                8.968703936700860024e-02,
                1.089329771368782868e-01,
                1.274352104115183404e-01,
                1.442232059743190353e-01,
                1.585209605389710785e-01,
                1.697992206073497823e-01,
                1.777863736909720827e-01,
                1.824469864020708543e-01,
                1.839397205857211670e-01,
                1.825658985584724392e-01,
                1.787176730860912666e-01,
                1.728317917757161704e-01,
                1.653521494452090401e-01,
                1.567021186831329960e-01,
                1.472661576201773348e-01,
                1.373793501883596724e-01,
                1.273231900217912382e-01,
                1.173259115464438673e-01,
                1.075658589701215528e-01,
                9.817665938293082539e-02,
                8.925325925656046877e-02,
                8.085815309826707709e-02,
                7.302735923472655766e-02,
                6.577587511923496666e-02,
                5.910247579657156719e-02,
                5.299381084138009290e-02,
                4.742781513856112602e-02,
                4.237648499479453046e-02,
                3.780808995871325767e-02,
                3.368889781496139935e-02,
                2.998448967748412047e-02,
                2.666073700370180655e-02,
            ];
            ddouble[] expected_dist_a3b4 = [
                6.291310563272764200e-07,
                1.524739335723873778e-06,
                3.489093675927508496e-06,
                7.564898120069278977e-06,
                1.559144296915498229e-05,
                3.064048680118335282e-05,
                5.758159646125797152e-05,
                1.037593184397517070e-04,
                1.797344278328408208e-04,
                3.000092224385862748e-04,
                4.836293022003270597e-04,
                7.545389669260120415e-04,
                1.141570355031978814e-03,
                1.677970616161966050e-03,
                2.400415380369625481e-03,
                3.347514921169945229e-03,
                4.557882911192471381e-03,
                6.067896427616657935e-03,
                7.909320551502688201e-03,
                1.010699477268236189e-02,
                1.267677840249518267e-02,
                1.562392919444086960e-02,
                1.894204766074140261e-02,
                2.261266573213621331e-02,
                2.660549997522472510e-02,
                3.087933376399200694e-02,
                3.538344547967899045e-02,
                4.005946494029297839e-02,
                4.484351968350430012e-02,
                4.966852634673985661e-02,
                5.446648856843914338e-02,
                5.917067891813388575e-02,
                6.371760520575917019e-02,
                6.804868795301644013e-02,
                7.211160298715951766e-02,
                7.586126874189065561e-02,
                7.926048026948553926e-02,
                8.228021017986784924e-02,
                8.489961030367489114e-02,
                8.710575687573439418e-02,
                8.889318684548604133e-02,
                9.026327415445103974e-02,
                9.122349320103542714e-02,
                9.178661298518313305e-02,
                9.196986029286058351e-02,
                9.179408436198104038e-02,
                9.128294927923621960e-02,
                9.046217428941430438e-02,
                8.935883654304563328e-02,
                8.800074575365095242e-02,
                8.641589588785808518e-02,
                8.463199540870405824e-02,
                8.267607472260452006e-02,
                8.057416729631071417e-02,
                7.835105934156649798e-02,
                7.603010192976084047e-02,
                7.363307881008866740e-02,
                7.118012297772959018e-02,
                6.868967509417983619e-02,
                6.617847712915700398e-02,
                6.366159501089561912e-02,
                6.115246458741532137e-02,
                5.866295577322193366e-02,
                5.620345035044795046e-02,
            ];

            foreach ((GumbelDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a1b2, expected_dist_a1b2), (dist_a2b2, expected_dist_a2b2),
                (dist_a3b4, expected_dist_a3b4),
            }) {
                for ((ddouble x, int i) = (-8, 0); i < expecteds.Length; x += 0.25, i++) {
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
            ddouble[] expected_dist_a1b1 = [
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
                1.357247607325042805e-289,
                1.071244712739173892e-225,
                6.210136486566144466e-176,
                3.537773075543882738e-137,
                5.390686197260364776e-107,
                1.730299058708983015e-83,
                3.507389196464635223e-65,
                6.344290125215141007e-51,
                8.054834089740903198e-40,
                3.577719320634459691e-31,
                1.942337604956407264e-24,
                3.414512624812977192e-19,
                4.150896920109045253e-15,
                6.300828916156514841e-12,
                1.892178694838292448e-09,
                1.608760113988778247e-07,
                5.119294298670732903e-06,
                7.577547728260715257e-05,
                6.179789893310934324e-04,
                3.168165149053242969e-03,
                1.131428638045962713e-02,
                3.049041346306221104e-02,
                6.598803584531254263e-02,
                1.203922620798295734e-01,
                1.922956455479649107e-01,
                2.769203340999089602e-01,
                3.678794411714423340e-01,
                4.589560693076638054e-01,
                5.452392118926050468e-01,
                6.235249162568003989e-01,
                6.922006275553463928e-01,
                7.508834766393948090e-01,
                8.000107130043535575e-01,
                8.404868737475783558e-01,
                8.734230184931166541e-01,
                8.999651626606277599e-01,
                9.211936551755157687e-01,
                9.380726685202481763e-01,
                9.514319929004534382e-01,
                9.619678894422100113e-01,
                9.702540025910624255e-01,
                9.767566411323356235e-01,
                9.818510730616665239e-01,
                9.858370182755026301e-01,
                9.889524805037951394e-01,
                9.913856230123334612e-01,
                9.932847020678414740e-01,
                9.947662257740502723e-01,
                9.959215680475387300e-01,
                9.968222788809181223e-01,
                9.975243173927524909e-01,
                9.980714079919146275e-01,
                9.984976904055608005e-01,
                9.988298055912923079e-01,
            ];
            ddouble[] expected_dist_a2b1 = [
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
                1.357247607325042805e-289,
                1.071244712739173892e-225,
                6.210136486566144466e-176,
                3.537773075543882738e-137,
                5.390686197260364776e-107,
                1.730299058708983015e-83,
                3.507389196464635223e-65,
                6.344290125215141007e-51,
                8.054834089740903198e-40,
                3.577719320634459691e-31,
                1.942337604956407264e-24,
                3.414512624812977192e-19,
                4.150896920109045253e-15,
                6.300828916156514841e-12,
                1.892178694838292448e-09,
                1.608760113988778247e-07,
                5.119294298670732903e-06,
                7.577547728260715257e-05,
                6.179789893310934324e-04,
                3.168165149053242969e-03,
                1.131428638045962713e-02,
                3.049041346306221104e-02,
                6.598803584531254263e-02,
                1.203922620798295734e-01,
                1.922956455479649107e-01,
                2.769203340999089602e-01,
                3.678794411714423340e-01,
                4.589560693076638054e-01,
                5.452392118926050468e-01,
                6.235249162568003989e-01,
                6.922006275553463928e-01,
                7.508834766393948090e-01,
                8.000107130043535575e-01,
                8.404868737475783558e-01,
                8.734230184931166541e-01,
                8.999651626606277599e-01,
                9.211936551755157687e-01,
                9.380726685202481763e-01,
                9.514319929004534382e-01,
                9.619678894422100113e-01,
                9.702540025910624255e-01,
                9.767566411323356235e-01,
                9.818510730616665239e-01,
                9.858370182755026301e-01,
                9.889524805037951394e-01,
                9.913856230123334612e-01,
                9.932847020678414740e-01,
                9.947662257740502723e-01,
                9.959215680475387300e-01,
                9.968222788809181223e-01,
            ];
            ddouble[] expected_dist_a1b2 = [
                8.054834089740903198e-40,
                3.160210699807487323e-35,
                3.577719320634459691e-31,
                1.352545348371412039e-27,
                1.942337604956407264e-24,
                1.187189396386732698e-21,
                3.414512624812977192e-19,
                5.049167717544199166e-17,
                4.150896920109045253e-15,
                2.032613549287901939e-13,
                6.300828916156514841e-12,
                1.304677438179973945e-10,
                1.892178694838292448e-09,
                2.004223336737291293e-08,
                1.608760113988778247e-07,
                1.010996392615728130e-06,
                5.119294298670732903e-06,
                2.142369113111383411e-05,
                7.577547728260715257e-05,
                2.310451324101896729e-04,
                6.179789893310934324e-04,
                1.472462484282561004e-03,
                3.168165149053242969e-03,
                6.229750238093417422e-03,
                1.131428638045962713e-02,
                1.915719869702440029e-02,
                3.049041346306221104e-02,
                4.594929151446972698e-02,
                6.598803584531254263e-02,
                9.082004171774372969e-02,
                1.203922620798295734e-01,
                1.543942385108051374e-01,
                1.922956455479649107e-01,
                2.334023677409892417e-01,
                2.769203340999089602e-01,
                3.220178007714426527e-01,
                3.678794411714423340e-01,
                4.137485310685763418e-01,
                4.589560693076638054e-01,
                5.029375464134954177e-01,
                5.452392118926050468e-01,
                5.855161995016551213e-01,
                6.235249162568003989e-01,
                6.591118581362204187e-01,
                6.922006275553463928e-01,
                7.227784991200458808e-01,
                7.508834766393948090e-01,
                7.765924405190635094e-01,
                8.000107130043535575e-01,
                8.212631680109220289e-01,
                8.404868737475783558e-01,
                8.578251682679061130e-01,
                8.734230184931166541e-01,
                8.874234913601654062e-01,
                8.999651626606277599e-01,
                9.111802979172057837e-01,
                9.211936551755157687e-01,
                9.301217782669365386e-01,
                9.380726685202481763e-01,
                9.451457415627220193e-01,
                9.514319929004534382e-01,
                9.570143109948293647e-01,
                9.619678894422100113e-01,
                9.663607006837475755e-01,
            ];
            ddouble[] expected_dist_a2b2 = [
                3.507389196464635223e-65,
                1.314159226736915273e-57,
                6.344290125215141007e-51,
                5.020633404967196589e-45,
                8.054834089740903198e-40,
                3.160210699807487323e-35,
                3.577719320634459691e-31,
                1.352545348371412039e-27,
                1.942337604956407264e-24,
                1.187189396386732698e-21,
                3.414512624812977192e-19,
                5.049167717544199166e-17,
                4.150896920109045253e-15,
                2.032613549287901939e-13,
                6.300828916156514841e-12,
                1.304677438179973945e-10,
                1.892178694838292448e-09,
                2.004223336737291293e-08,
                1.608760113988778247e-07,
                1.010996392615728130e-06,
                5.119294298670732903e-06,
                2.142369113111383411e-05,
                7.577547728260715257e-05,
                2.310451324101896729e-04,
                6.179789893310934324e-04,
                1.472462484282561004e-03,
                3.168165149053242969e-03,
                6.229750238093417422e-03,
                1.131428638045962713e-02,
                1.915719869702440029e-02,
                3.049041346306221104e-02,
                4.594929151446972698e-02,
                6.598803584531254263e-02,
                9.082004171774372969e-02,
                1.203922620798295734e-01,
                1.543942385108051374e-01,
                1.922956455479649107e-01,
                2.334023677409892417e-01,
                2.769203340999089602e-01,
                3.220178007714426527e-01,
                3.678794411714423340e-01,
                4.137485310685763418e-01,
                4.589560693076638054e-01,
                5.029375464134954177e-01,
                5.452392118926050468e-01,
                5.855161995016551213e-01,
                6.235249162568003989e-01,
                6.591118581362204187e-01,
                6.922006275553463928e-01,
                7.227784991200458808e-01,
                7.508834766393948090e-01,
                7.765924405190635094e-01,
                8.000107130043535575e-01,
                8.212631680109220289e-01,
                8.404868737475783558e-01,
                8.578251682679061130e-01,
                8.734230184931166541e-01,
                8.874234913601654062e-01,
                8.999651626606277599e-01,
                9.111802979172057837e-01,
                9.211936551755157687e-01,
                9.301217782669365386e-01,
                9.380726685202481763e-01,
                9.451457415627220193e-01,
            ];
            ddouble[] expected_dist_a3b4 = [
                1.608760113988778247e-07,
                4.150392558461454025e-07,
                1.010996392615728130e-06,
                2.333369221694511134e-06,
                5.119294298670732903e-06,
                1.070934359536413095e-05,
                2.142369113111383411e-05,
                4.109425507290971700e-05,
                7.577547728260715257e-05,
                1.346404075401420669e-04,
                2.310451324101896729e-04,
                3.837154583436779747e-04,
                6.179789893310934324e-04,
                9.669383473160519807e-04,
                1.472462484282561004e-03,
                2.185867155014639338e-03,
                3.168165149053242969e-03,
                4.489791902692240322e-03,
                6.229750238093417422e-03,
                8.474164776923377371e-03,
                1.131428638045962713e-02,
                1.484403383901214575e-02,
                1.915719869702440029e-02,
                2.434446566735199921e-02,
                3.049041346306221104e-02,
                3.767065882918071290e-02,
                4.594929151446972698e-02,
                5.537672254618709877e-02,
                6.598803584531254263e-02,
                7.780189754451265205e-02,
                9.082004171774372969e-02,
                1.050273184223726869e-01,
                1.203922620798295734e-01,
                1.368681165690696644e-01,
                1.543942385108051374e-01,
                1.728977919657045115e-01,
                1.922956455479649107e-01,
                2.124963858272356787e-01,
                2.334023677409892417e-01,
                2.549117324060547474e-01,
                2.769203340999089602e-01,
                2.993235303053610741e-01,
                3.220178007714426527e-01,
                3.449021729187441432e-01,
                3.678794411714423340e-01,
                3.908571766590693830e-01,
                4.137485310685763418e-01,
                4.364728442251377682e-01,
                4.589560693076638054e-01,
                4.811310325907158791e-01,
                5.029375464134954177e-01,
                5.243223948914627064e-01,
                5.452392118926050468e-01,
                5.656482701763533294e-01,
                5.855161995016551213e-01,
                6.048156500955050863e-01,
                6.235249162568003989e-01,
                6.416275331537114601e-01,
                6.591118581362204187e-01,
                6.759706461906954678e-01,
                6.922006275553463928e-01,
                7.078020940250089321e-01,
                7.227784991200458808e-01,
                7.371360760869097861e-01,
            ];

            foreach ((GumbelDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a1b2, expected_dist_a1b2), (dist_a2b2, expected_dist_a2b2),
                (dist_a3b4, expected_dist_a3b4),
            }) {
                for ((ddouble x, int i) = (-8, 0); i < expecteds.Length; x += 0.25, i++) {
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