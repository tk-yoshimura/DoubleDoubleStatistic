using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class NoncentralChiSquareDistributionTests {
        readonly NoncentralChiSquareDistribution dist_nu1lambda1 = new(nu: 1, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu2lambda1 = new(nu: 2, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu3lambda1 = new(nu: 3, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu4lambda1 = new(nu: 4, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu8lambda1 = new(nu: 8, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu16lambda1 = new(nu: 16, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu32lambda1 = new(nu: 32, lambda: 1);

        readonly NoncentralChiSquareDistribution dist_nu1lambda2 = new(nu: 1, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu2lambda2 = new(nu: 2, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu3lambda2 = new(nu: 3, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu4lambda2 = new(nu: 4, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu8lambda2 = new(nu: 8, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu16lambda2 = new(nu: 16, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu32lambda2 = new(nu: 32, lambda: 2);

        readonly NoncentralChiSquareDistribution dist_nu1lambda4 = new(nu: 1, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu2lambda4 = new(nu: 2, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu3lambda4 = new(nu: 3, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu4lambda4 = new(nu: 4, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu8lambda4 = new(nu: 8, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu16lambda4 = new(nu: 16, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu32lambda4 = new(nu: 32, lambda: 4);

        readonly NoncentralChiSquareDistribution dist_nu1lambda8 = new(nu: 1, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu2lambda8 = new(nu: 2, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu3lambda8 = new(nu: 3, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu4lambda8 = new(nu: 4, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu8lambda8 = new(nu: 8, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu16lambda8 = new(nu: 16, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu32lambda8 = new(nu: 32, lambda: 8);

        readonly NoncentralChiSquareDistribution dist_nu1lambda16 = new(nu: 1, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu2lambda16 = new(nu: 2, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu3lambda16 = new(nu: 3, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu4lambda16 = new(nu: 4, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu8lambda16 = new(nu: 8, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu16lambda16 = new(nu: 16, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu32lambda16 = new(nu: 32, lambda: 16);

        readonly NoncentralChiSquareDistribution dist_nu1lambda32 = new(nu: 1, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu2lambda32 = new(nu: 2, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu3lambda32 = new(nu: 3, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu4lambda32 = new(nu: 4, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu8lambda32 = new(nu: 8, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu16lambda32 = new(nu: 16, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu32lambda32 = new(nu: 32, lambda: 32);


        NoncentralChiSquareDistribution[] Dists => [
            dist_nu1lambda1,
            dist_nu2lambda1,
            dist_nu3lambda1,
            dist_nu4lambda1,
            dist_nu8lambda1,
            dist_nu16lambda1,
            dist_nu32lambda1,
            dist_nu1lambda2,
            dist_nu2lambda2,
            dist_nu3lambda2,
            dist_nu4lambda2,
            dist_nu8lambda2,
            dist_nu16lambda2,
            dist_nu32lambda2,
            dist_nu1lambda4,
            dist_nu2lambda4,
            dist_nu3lambda4,
            dist_nu4lambda4,
            dist_nu8lambda4,
            dist_nu16lambda4,
            dist_nu32lambda4,
            dist_nu1lambda8,
            dist_nu2lambda8,
            dist_nu3lambda8,
            dist_nu4lambda8,
            dist_nu8lambda8,
            dist_nu16lambda8,
            dist_nu32lambda8,
            dist_nu1lambda16,
            dist_nu2lambda16,
            dist_nu3lambda16,
            dist_nu4lambda16,
            dist_nu8lambda16,
            dist_nu16lambda16,
            dist_nu32lambda16,
            dist_nu1lambda32,
            dist_nu2lambda32,
            dist_nu3lambda32,
            dist_nu4lambda32,
            dist_nu8lambda32,
            dist_nu16lambda32,
            dist_nu32lambda32,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Console.WriteLine(dist.Mode);

                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode - 1e-10), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode + 1e-10), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536 * 2);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (NoncentralChiSquareDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 5) * 0.1, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            ddouble[] expected_nu1lambda1 = [
                ddouble.PositiveInfinity,
                6.835345838953019726e-01,
                4.815829224301914624e-01,
                3.909444553678056433e-01,
                3.359531306989162469e-01,
                2.976251561624582154e-01,
                2.686576484629918382e-01,
                2.455695740489495238e-01,
                2.264666234573106118e-01,
                2.102198165046287870e-01,
                1.961098871858855952e-01,
                1.836544906192093818e-01,
                1.725166836528668202e-01,
                1.624530889920689503e-01,
                1.532829138795417656e-01,
                1.448685981032288295e-01,
                1.371032722715034569e-01,
                1.299023716065130785e-01,
                1.231978757526228097e-01,
                1.169342589457169901e-01,
                1.110655838723473426e-01,
                1.055533781928278819e-01,
                1.003650577536092381e-01,
                9.547273870139739993e-02,
                9.085233082365826163e-02,
                8.648283727993190395e-02,
                8.234580784466316583e-02,
                7.842490773097611501e-02,
                7.470557441219075057e-02,
                7.117474212931819422e-02,
                6.782061895401080009e-02,
                6.463250501575519558e-02,
                6.160064323277038983e-02,
                5.871609590206318396e-02,
                5.597064200760477198e-02,
                5.335669123752823095e-02,
                5.086721156095376506e-02,
                4.849566787338715351e-02,
                4.623596972766299412e-02,
                4.408242656226735862e-02,
                4.202970914787443335e-02,
                4.007281621627391960e-02,
                3.820704542866273229e-02,
                3.642796799388201001e-02,
                3.473140637022924943e-02,
                3.311341458353388167e-02,
                3.157026077432105449e-02,
                3.009841165200597204e-02,
                2.869451858720413126e-02,
                2.735540511679312656e-02,
                2.607805567218754797e-02,
                2.485960537086998393e-02,
                2.369733073573474999e-02,
                2.258864122718264034e-02,
                2.153107148990614228e-02,
                2.052227423053128857e-02,
                1.956001365422160507e-02,
                1.864215939840018280e-02,
                1.776668091022762935e-02,
                1.693164222165328275e-02,
                1.613519708194883703e-02,
                1.537558441281628059e-02,
                1.465112405558203351e-02,
                1.396021278376823917e-02,
            ];
            ddouble[] expected_nu2lambda4 = [
                6.766764161830635116e-02,
                7.176563760889131949e-02,
                7.560500290056605677e-02,
                7.919170383617120390e-02,
                8.253196319794596503e-02,
                8.563222063419634322e-02,
                8.849909592719816531e-02,
                9.113935495467327674e-02,
                9.355987820265798671e-02,
                9.576763169299580547e-02,
                9.776964019401011841e-02,
                9.957296258814543610e-02,
                1.011846692754947707e-01,
                1.026118214971528941e-01,
                1.038614524672426959e-01,
                1.049405502072503771e-01,
                1.058560419809717601e-01,
                1.066147802329105998e-01,
                1.072235299373853062e-01,
                1.076889572698781911e-01,
                1.080176195163177200e-01,
                1.082159561400005271e-01,
                1.082902809297507790e-01,
                1.082467751566710268e-01,
                1.080914816704661441e-01,
                1.078302998698204007e-01,
                1.074689814846764013e-01,
                1.070131271115108418e-01,
                1.064681834458199167e-01,
                1.058394411590266176e-01,
                1.051320333698993748e-01,
                1.043509346633337276e-01,
                1.035009606119933201e-01,
                1.025867677588442467e-01,
                1.016128540210380488e-01,
                1.005835594779233150e-01,
                9.950306750817747292e-02,
                9.837540624316981619e-02,
                9.720445030568361511e-02,
                9.599392280505165964e-02,
                9.474739756159156201e-02,
                9.346830153497398541e-02,
                9.215991743281641557e-02,
                9.082538647737445348e-02,
                8.946771130970025687e-02,
                8.808975901206167469e-02,
                8.669426423076474764e-02,
                8.528383238280086431e-02,
                8.386094293095089358e-02,
                8.242795271312433614e-02,
                8.098709931279596264e-02,
                7.954050445842328343e-02,
                7.809017744069841893e-02,
                7.663801853739518610e-02,
                7.518582243643325003e-02,
                7.373528164858955680e-02,
                7.228798990204740593e-02,
                7.084544551168871218e-02,
                6.940905471670615057e-02,
                6.798013498074247096e-02,
                6.655991824935311718e-02,
                6.514955416014200917e-02,
                6.375011320143617044e-02,
                6.236258981584837102e-02,
            ];
            ddouble[] expected_nu4lambda2 = [
                0.000000000000000000e+00,
                1.114073588566949491e-02,
                2.158565013451052908e-02,
                3.135773439302179388e-02,
                4.048019857076348010e-02,
                4.897636319596097160e-02,
                5.686956274970249303e-02,
                6.418305915097749503e-02,
                7.093996461786045149e-02,
                7.716317318020474159e-02,
                8.287530016641460329e-02,
                8.809862903136270185e-02,
                9.285506493444331544e-02,
                9.716609451623721416e-02,
                1.010527513594318832e-01,
                1.045355866546173601e-01,
                1.076346446244688115e-01,
                1.103694422907527811e-01,
                1.127589531976506410e-01,
                1.148215947321828551e-01,
                1.165752187081383756e-01,
                1.180371049039534298e-01,
                1.192239572675274362e-01,
                1.201519025221000903e-01,
                1.208364909271112675e-01,
                1.212926989665097793e-01,
                1.215349337543261271e-01,
                1.215770389635619370e-01,
                1.214323020996189162e-01,
                1.211134629536729096e-01,
                1.206327230846400839e-01,
                1.200017561907346791e-01,
                1.192317192431484874e-01,
                1.183332642651208771e-01,
                1.173165506496740007e-01,
                1.161912579186043237e-01,
                1.149665988339816092e-01,
                1.136513327814619306e-01,
                1.122537793521990734e-01,
                1.107818320570792270e-01,
                1.092429721134428600e-01,
                1.076442822504171531e-01,
                1.059924604845072782e-01,
                1.042938338221947175e-01,
                1.025543718510127933e-01,
                1.007797001849175256e-01,
                9.897511373378979627e-02,
                9.714558977059700728e-02,
                9.529580077314514508e-02,
                9.343012702047036322e-02,
                9.155266892678894097e-02,
                8.966725909854246501e-02,
                8.777747410247942195e-02,
                8.588664593490157040e-02,
                8.399787318420250115e-02,
                8.211403188063969016e-02,
                8.023778602893161249e-02,
                7.837159782076301429e-02,
                7.651773752562934794e-02,
                7.467829305966718068e-02,
                7.285517923320702915e-02,
                7.105014667876045908e-02,
                6.926479046202510537e-02,
                6.750055837926194346e-02,
            ];

            foreach ((NoncentralChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1lambda1, expected_nu1lambda1),
                (dist_nu2lambda4, expected_nu2lambda4),
                (dist_nu4lambda2, expected_nu4lambda2),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else if (ddouble.IsPositiveInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_nu1lambda1 = [
                0.000000000000000000e+00,
                1.710556289337611113e-01,
                2.417303374571271590e-01,
                2.957053740472119219e-01,
                3.409007826677720732e-01,
                3.803748500113006692e-01,
                4.156924256045555488e-01,
                4.477822188423997152e-01,
                4.772498680518171277e-01,
                5.045173382944607710e-01,
                5.298935418727288571e-01,
                5.536135530746811861e-01,
                5.758619373747236869e-01,
                5.967873977330609980e-01,
                6.165123764677900020e-01,
                6.351395783764266501e-01,
                6.527565366822647430e-01,
                6.694388914556553827e-01,
                6.852527959482315190e-01,
                7.002567171631854936e-01,
                7.145028063106916294e-01,
                7.280379579077228547e-01,
                7.409046395947428421e-01,
                7.531415505018896894e-01,
                7.647841496310278098e-01,
                7.758650844536137958e-01,
                7.864145420317905089e-01,
                7.964605393527692812e-01,
                8.060291655107413300e-01,
                8.151447854031753648e-01,
                8.238302124111297475e-01,
                8.321068558880315713e-01,
                8.399948480369093806e-01,
                8.475131538057446567e-01,
                8.546796666985161384e-01,
                8.615112928309222040e-01,
                8.680240251147774977e-01,
                8.742330091043546414e-01,
                8.801526017596598139e-01,
                8.857964241594525889e-01,
                8.911774090183970376e-01,
                8.963078437186245973e-01,
                9.011994094490448637e-01,
                9.058632169502979226e-01,
                9.103098392850577092e-01,
                9.145493419889939579e-01,
                9.185913109044827607e-01,
                9.224448779549541122e-01,
                9.261187450809619248e-01,
                9.296212065282319914e-01,
                9.329601696521022669e-01,
                9.361431743809865980e-01,
                9.391774114630597747e-01,
                9.420697396047645888e-01,
                9.448267015964423310e-01,
                9.474545395090433475e-01,
                9.499592090361757002e-01,
                9.523463930473931027e-01,
                9.546215144114523765e-01,
                9.567897481420460926e-01,
                9.588560329131549098e-01,
                9.608250819864625480e-01,
                9.627013935892110297e-01,
                9.644892607773055682e-01,
            ];
            ddouble[] expected_nu2lambda4 = [
                0.000000000000000000e+00,
                8.717304152254530536e-03,
                1.793063270833497938e-02,
                2.760802684022219122e-02,
                3.771828996773488041e-02,
                4.823101730973757689e-02,
                5.911662066782657349e-02,
                7.034634878836100447e-02,
                8.189230363059368800e-02,
                9.372745285122097147e-02,
                1.058256387989098407e-01,
                1.181615842962573659e-01,
                1.307108954711830862e-01,
                1.434500618849348652e-01,
                1.563564541896880211e-01,
                1.694083195351262505e-01,
                1.825847749303874668e-01,
                1.958657987553330537e-01,
                2.092322206032292409e-01,
                2.226657096256068058e-01,
                2.361487615392676709e-01,
                2.496646844451191893e-01,
                2.631975835987281909e-01,
                2.767323452631686997e-01,
                2.902546197658875116e-01,
                3.037508038729033744e-01,
                3.172080225856778735e-01,
                3.306141104584294221e-01,
                3.439575925264986567e-01,
                3.572276649295851314e-01,
                3.704141753072602294e-01,
                3.835076030380932810e-01,
                3.964990393880056252e-01,
                4.093801676280598434e-01,
                4.221432431767996318e-01,
                4.347810738174570644e-01,
                4.472870000358389331e-01,
                4.596548755204474257e-01,
                4.718790478624179108e-01,
                4.839543394891167027e-01,
                4.958760288617327650e-01,
                5.076398319639260359e-01,
                5.192418841055250756e-01,
                5.306787220624030743e-01,
                5.419472665709849091e-01,
                5.530448051933655185e-01,
                5.639689755666781323e-01,
                5.747177490482244799e-01,
                5.852894147658687096e-01,
                5.956825640813490086e-01,
                6.058960754724582731e-01,
                6.159290998384607807e-01,
                6.257810462316640843e-01,
                6.354515680167206471e-01,
                6.449405494580200404e-01,
                6.542480927344028752e-01,
                6.633745053794017910e-01,
                6.723202881442987833e-01,
                6.810861232804235055e-01,
                6.896728632363683387e-01,
                6.980815197651076831e-01,
                7.063132534353846470e-01,
                7.143693635411902942e-01,
                7.222512784026695609e-01,
            ];
            ddouble[] expected_nu4lambda2 = [
                0.000000000000000000e+00,
                7.036628899837029175e-04,
                2.756190936206948713e-03,
                6.072040539497646112e-03,
                1.056855687927373694e-02,
                1.616599434979718031e-02,
                2.278752423594962440e-02,
                3.035923094855051338e-02,
                3.881009803911745354e-02,
                4.807198512029116355e-02,
                5.807959673054680083e-02,
                6.877044409994319263e-02,
                8.008480069719130234e-02,
                9.196565236694045908e-02,
                1.043586427995838867e-01,
                1.172120150138039746e-01,
                1.304765494742261556e-01,
                1.441054994126964306e-01,
                1.580545238715859446e-01,
                1.722816189409333842e-01,
                1.867470476179581285e-01,
                2.014132686772898140e-01,
                2.162448649029930781e-01,
                2.312084709989254971e-01,
                2.462727014619801591e-01,
                2.614080786731791117e-01,
                2.765869614342371663e-01,
                2.917834741519695352e-01,
                3.069734368496291332e-01,
                3.221342961628153589e-01,
                3.372450574578507165e-01,
                3.522862181923884228e-01,
                3.672397026213696791e-01,
                3.820887979361884401e-01,
                3.968180919109701166e-01,
                4.114134121171114655e-01,
                4.258617667556117148e-01,
                4.401512871461251142e-01,
                4.542711719020692818e-01,
                4.682116328123970272e-01,
                4.819638424427683532e-01,
                4.955198834617606796e-01,
                5.088726996913587319e-01,
                5.220160488752614292e-01,
                5.349444571534163195e-01,
                5.476531752266708963e-01,
                5.601381361913944090e-01,
                5.723959150204013246e-01,
                5.844236896634041889e-01,
                5.962192037375420028e-01,
                6.077807307762206346e-01,
                6.191070400025268095e-01,
                6.301973635918273198e-01,
                6.410513653867927486e-01,
                6.516691110269787757e-01,
                6.620510394542344912e-01,
                6.721979357545417866e-01,
                6.821109052964531294e-01,
                6.917913491259911618e-01,
                7.012409405777811910e-01,
                7.104616030621927480e-01,
                7.194554889884299254e-01,
                7.282249597837771038e-01,
                7.367725669695810442e-01,
            ];

            foreach ((NoncentralChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1lambda1, expected_nu1lambda1),
                (dist_nu2lambda4, expected_nu2lambda4),
                (dist_nu4lambda2, expected_nu4lambda2),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }
    }
}