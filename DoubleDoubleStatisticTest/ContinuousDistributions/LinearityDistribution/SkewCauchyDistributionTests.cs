using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class SkewCauchyDistributionTests {
        readonly SkewCauchyDistribution dist_alpha000mu1gamma1 = new(alpha: 0, mu: 1, gamma: 1);
        readonly SkewCauchyDistribution dist_alpha025mu1gamma1 = new(alpha: 0.25, mu: 1, gamma: 1);
        readonly SkewCauchyDistribution dist_alpha050mu1gamma2 = new(alpha: 0.50, mu: 1, gamma: 2);
        readonly SkewCauchyDistribution dist_alpha025mu2gamma1 = new(alpha: 0.25, mu: 2, gamma: 1);
        readonly SkewCauchyDistribution dist_alpha050mu2gamma2 = new(alpha: 0.50, mu: 2, gamma: 2);
        readonly SkewCauchyDistribution dist_alpha075mu4gamma3 = new(alpha: 0.75, mu: 4, gamma: 3);
        readonly SkewCauchyDistribution dist_alpham050mu1gamma2 = new(alpha: -0.50, mu: 1, gamma: 2);

        SkewCauchyDistribution[] Dists => [
            dist_alpha000mu1gamma1,
            dist_alpha025mu1gamma1,
            dist_alpha050mu1gamma2,
            dist_alpha025mu2gamma1,
            dist_alpha050mu2gamma2,
            dist_alpha075mu4gamma3,
            dist_alpham050mu1gamma2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Gamma={dist.Gamma}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
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

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-280"; p /= 10) {
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
            foreach (SkewCauchyDistribution dist in Dists) {
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
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (SkewCauchyDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 2) * 0.05, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (SkewCauchyDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                (SkewCauchyDistribution? dist_fit, ddouble error) = SkewCauchyDistribution.Fit(xs, (0.2, 0.8));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (SkewCauchyDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha000mu1gamma1 = [
                1.224268793014579408e-02,
                1.350917288843674069e-02,
                1.497928876159014891e-02,
                1.669822353751032953e-02,
                1.872411095198768527e-02,
                2.113260655162095547e-02,
                2.402338763651250264e-02,
                2.752950366994946774e-02,
                3.183098861837906773e-02,
                3.717487721854489252e-02,
                4.390481188741940377e-02,
                5.250472349423351837e-02,
                6.366197723675813547e-02,
                7.835320275293310155e-02,
                9.794150344116635265e-02,
                1.242184921692841837e-01,
                1.591549430918953456e-01,
                2.037183271576260279e-01,
                2.546479089470325419e-01,
                2.995857752318029643e-01,
                3.183098861837906912e-01,
                2.995857752318029643e-01,
                2.546479089470325419e-01,
                2.037183271576260279e-01,
                1.591549430918953456e-01,
                1.242184921692841837e-01,
                9.794150344116635265e-02,
                7.835320275293310155e-02,
                6.366197723675813547e-02,
                5.250472349423351837e-02,
                4.390481188741940377e-02,
                3.717487721854489252e-02,
                3.183098861837906773e-02,
            ];
            ddouble[] expected_dist_alpha025mu1gamma1 = [
                7.004374023604195698e-03,
                7.742672907173286066e-03,
                8.602969896859206933e-03,
                9.613385824342670846e-03,
                1.081052443643062619e-02,
                1.224268793014579408e-02,
                1.397458036904446893e-02,
                1.609432008794447339e-02,
                1.872411095198768527e-02,
                2.203683827426243177e-02,
                2.628246766655152611e-02,
                3.183098861837906773e-02,
                3.924368459800159359e-02,
                4.939291337334682230e-02,
                6.366197723675813547e-02,
                8.425849928394460453e-02,
                1.145915590261646633e-01,
                1.591549430918953456e-01,
                2.203683827426242969e-01,
                2.864788975654116165e-01,
                3.183098861837906912e-01,
                3.060671982536448676e-01,
                2.744050742963712874e-01,
                2.340513868998461144e-01,
                1.940913940145065075e-01,
                1.591549430918953456e-01,
                1.304548713867994536e-01,
                1.075371237107400901e-01,
                8.941288937746928733e-02,
                7.507308636410156555e-02,
                6.366197723675813547e-02,
                5.450511749722443400e-02,
                4.708726126979152077e-02,
            ];
            ddouble[] expected_dist_alpha050mu1gamma2 = [
                6.121343965072897041e-03,
                6.754586444218370345e-03,
                7.489644380795074455e-03,
                8.349111768755164767e-03,
                9.362055475993842635e-03,
                1.056630327581047774e-02,
                1.201169381825625132e-02,
                1.376475183497473387e-02,
                1.591549430918953387e-02,
                1.858743860927244626e-02,
                2.195240594370970189e-02,
                2.625236174711675918e-02,
                3.183098861837906773e-02,
                3.917660137646655077e-02,
                4.897075172058317633e-02,
                6.210924608464209185e-02,
                7.957747154594767280e-02,
                1.018591635788130140e-01,
                1.273239544735162709e-01,
                1.497928876159014822e-01,
                1.591549430918953456e-01,
                1.580573227947098536e-01,
                1.548534581434657387e-01,
                1.497928876159014822e-01,
                1.432394487827058083e-01,
                1.356113124569995609e-01,
                1.273239544735162709e-01,
                1.187477295607923811e-01,
                1.101841913713121485e-01,
                1.018591635788130140e-01,
                9.392750739849560493e-02,
                8.648419549144500951e-02,
                7.957747154594767280e-02,
            ];
            ddouble[] expected_dist_alpha025mu2gamma1 = [
                4.897075172058318847e-03,
                5.324886571847799177e-03,
                5.810930985099627480e-03,
                6.366197723675813373e-03,
                7.004374023604195698e-03,
                7.742672907173286066e-03,
                8.602969896859206933e-03,
                9.613385824342670846e-03,
                1.081052443643062619e-02,
                1.224268793014579408e-02,
                1.397458036904446893e-02,
                1.609432008794447339e-02,
                1.872411095198768527e-02,
                2.203683827426243177e-02,
                2.628246766655152611e-02,
                3.183098861837906773e-02,
                3.924368459800159359e-02,
                4.939291337334682230e-02,
                6.366197723675813547e-02,
                8.425849928394460453e-02,
                1.145915590261646633e-01,
                1.591549430918953456e-01,
                2.203683827426242969e-01,
                2.864788975654116165e-01,
                3.183098861837906912e-01,
                3.060671982536448676e-01,
                2.744050742963712874e-01,
                2.340513868998461144e-01,
                1.940913940145065075e-01,
                1.591549430918953456e-01,
                1.304548713867994536e-01,
                1.075371237107400901e-01,
                8.941288937746928733e-02,
            ];
            ddouble[] expected_dist_alpha050mu2gamma2 = [
                4.301484948429603467e-03,
                4.672438696275826170e-03,
                5.092958178940651219e-03,
                5.572164309563074937e-03,
                6.121343965072897041e-03,
                6.754586444218370345e-03,
                7.489644380795074455e-03,
                8.349111768755164767e-03,
                9.362055475993842635e-03,
                1.056630327581047774e-02,
                1.201169381825625132e-02,
                1.376475183497473387e-02,
                1.591549430918953387e-02,
                1.858743860927244626e-02,
                2.195240594370970189e-02,
                2.625236174711675918e-02,
                3.183098861837906773e-02,
                3.917660137646655077e-02,
                4.897075172058317633e-02,
                6.210924608464209185e-02,
                7.957747154594767280e-02,
                1.018591635788130140e-01,
                1.273239544735162709e-01,
                1.497928876159014822e-01,
                1.591549430918953456e-01,
                1.580573227947098536e-01,
                1.548534581434657387e-01,
                1.497928876159014822e-01,
                1.432394487827058083e-01,
                1.356113124569995609e-01,
                1.273239544735162709e-01,
                1.187477295607923811e-01,
                1.101841913713121485e-01,
            ];
            ddouble[] expected_dist_alpha075mu4gamma3 = [
                9.244236772036515089e-04,
                9.844635655168784260e-04,
                1.050527677174226734e-03,
                1.123446657119261515e-03,
                1.204198812801225487e-03,
                1.293942626763376841e-03,
                1.394057895695433600e-03,
                1.506198199607842016e-03,
                1.632358390686106282e-03,
                1.774962190615932698e-03,
                1.936976995033209449e-03,
                2.122065907891937647e-03,
                2.334791341201398421e-03,
                2.580890969057762600e-03,
                2.867656632286402166e-03,
                3.204461941447557093e-03,
                3.603508145476875396e-03,
                4.080895976715264405e-03,
                4.658193456348156021e-03,
                5.364773362648157795e-03,
                6.241370317329228423e-03,
                7.345612758087476969e-03,
                8.760822555517173635e-03,
                1.061032953945968867e-02,
                1.308122819933386395e-02,
                1.646430445778227294e-02,
                2.122065907891937733e-02,
                2.808616642798153137e-02,
                3.819718634205488544e-02,
                5.305164769729844854e-02,
                7.345612758087476102e-02,
                9.549296585513721014e-02,
                1.061032953945968971e-01,

            ];
            ddouble[] expected_dist_alpham050mu1gamma2 = [
                4.212924964197230226e-02,
                4.538279565392659004e-02,
                4.897075172058317633e-02,
                5.292912657097673113e-02,
                5.729577951308233164e-02,
                6.210924608464209185e-02,
                6.740679942715567530e-02,
                7.322144346719786090e-02,
                7.957747154594767280e-02,
                8.648419549144500951e-02,
                9.392750739849560493e-02,
                1.018591635788130140e-01,
                1.101841913713121485e-01,
                1.187477295607923811e-01,
                1.273239544735162709e-01,
                1.356113124569995609e-01,
                1.432394487827058083e-01,
                1.497928876159014822e-01,
                1.548534581434657387e-01,
                1.580573227947098536e-01,
                1.591549430918953456e-01,
                1.497928876159014822e-01,
                1.273239544735162709e-01,
                1.018591635788130140e-01,
                7.957747154594767280e-02,
                6.210924608464209185e-02,
                4.897075172058317633e-02,
                3.917660137646655077e-02,
                3.183098861837906773e-02,
                2.625236174711675918e-02,
                2.195240594370970189e-02,
                1.858743860927244626e-02,
                1.591549430918953387e-02,
            ];

            foreach ((SkewCauchyDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha000mu1gamma1 , expected_dist_alpha000mu1gamma1 ),
                (dist_alpha025mu1gamma1 , expected_dist_alpha025mu1gamma1 ),
                (dist_alpha050mu1gamma2 , expected_dist_alpha050mu1gamma2 ),
                (dist_alpha025mu2gamma1 , expected_dist_alpha025mu2gamma1 ),
                (dist_alpha050mu2gamma2 , expected_dist_alpha050mu2gamma2 ),
                (dist_alpha075mu4gamma3 , expected_dist_alpha075mu4gamma3 ),
                (dist_alpham050mu1gamma2, expected_dist_alpham050mu1gamma2),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.25, i++) {
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
            ddouble[] expected_dist_alpha000mu1gamma1 = [
                6.283295818900114238e-02,
                6.604810022015533688e-02,
                6.960448727306395300e-02,
                7.355844397326222195e-02,
                7.797913037736925457e-02,
                8.295231765631966292e-02,
                8.858553278290470789e-02,
                9.501516093917983241e-02,
                1.024163823495666814e-01,
                1.110172584549998231e-01,
                1.211189415908434097e-01,
                1.331249387476565182e-01,
                1.475836176504332631e-01,
                1.652493405385678793e-01,
                1.871670418109988021e-01,
                2.147767125227227347e-01,
                2.500000000000000000e-01,
                2.951672353008665262e-01,
                3.524163823495667369e-01,
                4.220208696226306899e-01,
                5.000000000000000000e-01,
                5.779791303773693656e-01,
                6.475836176504332631e-01,
                7.048327646991334738e-01,
                7.500000000000000000e-01,
                7.852232874772773208e-01,
                8.128329581890012534e-01,
                8.347506594614321207e-01,
                8.524163823495667369e-01,
                8.668750612523434818e-01,
                8.788810584091566458e-01,
                8.889827415450002324e-01,
                8.975836176504332631e-01,
            ];
            ddouble[] expected_dist_alpha025mu1gamma1 = [
                3.554485670811718956e-02,
                3.738594422873497969e-02,
                3.942634253344007300e-02,
                4.169991583933890800e-02,
                4.424856365064638863e-02,
                4.712471864175088454e-02,
                5.039482115421706210e-02,
                5.414423663298545408e-02,
                5.848434778302696868e-02,
                6.356299459607406277e-02,
                6.958018430830675838e-02,
                7.681228676217505269e-02,
                8.565018841493110546e-02,
                9.666079380686742883e-02,
                1.106877132378250028e-01,
                1.290156522169730036e-01,
                1.536245735243501054e-01,
                1.875000000000000000e-01,
                2.346247186417509123e-01,
                2.981877132378249473e-01,
                3.750000000000000000e-01,
                4.535411977362515112e-01,
                5.263986769885542483e-01,
                5.900260870282884040e-01,
                6.434708906534034600e-01,
                6.875000000000000000e-01,
                7.235724229703806021e-01,
                7.532105708890668438e-01,
                7.777403943883089976e-01,
                7.982319159786310792e-01,
                8.155204779369584767e-01,
                8.302503140248851388e-01,
                8.429176045274970885e-01,
            ];
            ddouble[] expected_dist_alpha050mu1gamma2 = [
                3.141647909450057119e-02,
                3.302405011007766844e-02,
                3.480224363653197650e-02,
                3.677922198663111097e-02,
                3.898956518868462728e-02,
                4.147615882815983146e-02,
                4.429276639145235395e-02,
                4.750758046958991621e-02,
                5.120819117478334070e-02,
                5.550862922749991157e-02,
                6.055947079542170486e-02,
                6.656246937382825912e-02,
                7.379180882521663154e-02,
                8.262467026928393965e-02,
                9.358352090549940105e-02,
                1.073883562613613674e-01,
                1.250000000000000000e-01,
                1.475836176504332631e-01,
                1.762081911747833685e-01,
                2.110104348113153450e-01,
                2.500000000000000000e-01,
                2.896970140893848078e-01,
                3.288526850668801460e-01,
                3.669686955660539929e-01,
                4.036245735243501054e-01,
                4.384988745670035382e-01,
                4.713754264756498946e-01,
                5.021369763627439031e-01,
                5.307505627164981199e-01,
                5.572491470487002108e-01,
                5.817130924355433663e-01,
                6.042537256500071230e-01,
                6.250000000000000000e-01,
            ];
            ddouble[] expected_dist_alpha025mu2gamma1 = [
                2.968756812042416815e-02,
                3.096419987988541234e-02,
                3.235485841010554697e-02,
                3.387542647564989462e-02,
                3.554485670811718956e-02,
                3.738594422873497969e-02,
                3.942634253344007300e-02,
                4.169991583933890800e-02,
                4.424856365064638863e-02,
                4.712471864175088454e-02,
                5.039482115421706210e-02,
                5.414423663298545408e-02,
                5.848434778302696868e-02,
                6.356299459607406277e-02,
                6.958018430830675838e-02,
                7.681228676217505269e-02,
                8.565018841493110546e-02,
                9.666079380686742883e-02,
                1.106877132378250028e-01,
                1.290156522169730036e-01,
                1.536245735243501054e-01,
                1.875000000000000000e-01,
                2.346247186417509123e-01,
                2.981877132378249473e-01,
                3.750000000000000000e-01,
                4.535411977362515112e-01,
                5.263986769885542483e-01,
                5.900260870282884040e-01,
                6.434708906534034600e-01,
                6.875000000000000000e-01,
                7.235724229703806021e-01,
                7.532105708890668438e-01,
                7.777403943883089976e-01,
            ];
            ddouble[] expected_dist_alpha050mu2gamma2 = [
                2.628422835562668758e-02,
                2.740501928634545004e-02,
                2.862457352435004987e-02,
                2.995638296545163892e-02,
                3.141647909450057119e-02,
                3.302405011007766844e-02,
                3.480224363653197650e-02,
                3.677922198663111097e-02,
                3.898956518868462728e-02,
                4.147615882815983146e-02,
                4.429276639145235395e-02,
                4.750758046958991621e-02,
                5.120819117478334070e-02,
                5.550862922749991157e-02,
                6.055947079542170486e-02,
                6.656246937382825912e-02,
                7.379180882521663154e-02,
                8.262467026928393965e-02,
                9.358352090549940105e-02,
                1.073883562613613674e-01,
                1.250000000000000000e-01,
                1.475836176504332631e-01,
                1.762081911747833685e-01,
                2.110104348113153450e-01,
                2.500000000000000000e-01,
                2.896970140893848078e-01,
                3.288526850668801460e-01,
                3.669686955660539929e-01,
                4.036245735243501054e-01,
                4.384988745670035382e-01,
                4.713754264756498946e-01,
                5.021369763627439031e-01,
                5.307505627164981199e-01,
            ];
            ddouble[] expected_dist_alpha075mu4gamma3 = [
                7.438645892854423014e-03,
                7.677139099522459209e-03,
                7.931379357638379535e-03,
                8.202973769125695980e-03,
                8.493754953174170730e-03,
                8.805821869319316320e-03,
                9.141589798858332228e-03,
                9.503851961987400121e-03,
                9.895856040141384757e-03,
                1.032139995996179949e-02,
                1.078495280336851103e-02,
                1.129180882521663154e-02,
                1.184828556937239652e-02,
                1.246198140957831269e-02,
                1.314211417781334379e-02,
                1.389997194644629341e-02,
                1.474952121688212492e-02,
                1.570823954725028559e-02,
                1.679827371807235403e-02,
                1.804807887766181340e-02,
                1.949478259434231364e-02,
                2.118766486535801630e-02,
                2.319339476943557687e-02,
                2.560409558739167035e-02,
                2.855006280497703053e-02,
                3.222026460228914757e-02,
                3.689590441260831577e-02,
                4.300521740565765860e-02,
                5.120819117478336846e-02,
                6.250000000000000000e-02,
                7.820823954725029947e-02,
                9.939590441260831577e-02,
                1.250000000000000000e-01,
            ];
            ddouble[] expected_dist_alpham050mu1gamma2 = [
                2.580313044339460071e-01,
                2.689637026214802273e-01,
                2.807505627164982309e-01,
                2.934799414016059371e-01,
                3.072491470487002108e-01,
                3.221650687840840743e-01,
                3.383441220417038875e-01,
                3.559115829780122864e-01,
                3.750000000000000000e-01,
                3.957462743499929325e-01,
                4.182869075644566892e-01,
                4.427508529512998448e-01,
                4.692494372835018246e-01,
                4.978630236372560969e-01,
                5.286245735243501054e-01,
                5.615011254329964618e-01,
                5.963754264756498946e-01,
                6.330313044339460626e-01,
                6.711473149331198540e-01,
                7.103029859106151367e-01,
                7.500000000000000000e-01,
                7.889895651886846828e-01,
                8.237918088252166315e-01,
                8.524163823495667369e-01,
                8.750000000000000000e-01,
                8.926116437386386604e-01,
                9.064164790945006267e-01,
                9.173753297307161159e-01,
                9.262081911747833685e-01,
                9.334375306261717409e-01,
                9.394405292045783229e-01,
                9.444913707725001162e-01,
                9.487918088252166315e-01,
            ];

            foreach ((SkewCauchyDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha000mu1gamma1 , expected_dist_alpha000mu1gamma1 ),
                (dist_alpha025mu1gamma1 , expected_dist_alpha025mu1gamma1 ),
                (dist_alpha050mu1gamma2 , expected_dist_alpha050mu1gamma2 ),
                (dist_alpha025mu2gamma1 , expected_dist_alpha025mu2gamma1 ),
                (dist_alpha050mu2gamma2 , expected_dist_alpha050mu2gamma2 ),
                (dist_alpha075mu4gamma3 , expected_dist_alpha075mu4gamma3 ),
                (dist_alpham050mu1gamma2, expected_dist_alpham050mu1gamma2),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.25, i++) {
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