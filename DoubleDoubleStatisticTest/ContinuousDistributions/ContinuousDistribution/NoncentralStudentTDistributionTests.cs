using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class NoncentralStudentTDistributionTests {
        readonly NoncentralStudentTDistribution dist_nu1mu0 = new(nu: 1, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu2mu0 = new(nu: 2, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu3mu0 = new(nu: 3, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu4mu0 = new(nu: 4, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu5mu0 = new(nu: 5, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu8mu0 = new(nu: 8, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu16mu0 = new(nu: 16, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu32mu0 = new(nu: 32, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu64mu0 = new(nu: 64, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu128mu0 = new(nu: 128, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu129mu0 = new(nu: 129, mu: 0);

        readonly NoncentralStudentTDistribution dist_nu1mu1 = new(nu: 1, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu2mu1 = new(nu: 2, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu3mu1 = new(nu: 3, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu4mu1 = new(nu: 4, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu5mu1 = new(nu: 5, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu8mu1 = new(nu: 8, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu16mu1 = new(nu: 16, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu32mu1 = new(nu: 32, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu64mu1 = new(nu: 64, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu128mu1 = new(nu: 128, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu129mu1 = new(nu: 129, mu: 1);

        readonly NoncentralStudentTDistribution dist_nu1mum1 = new(nu: 1, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu2mum1 = new(nu: 2, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu3mum1 = new(nu: 3, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu4mum1 = new(nu: 4, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu5mum1 = new(nu: 5, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu8mum1 = new(nu: 8, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu16mum1 = new(nu: 16, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu32mum1 = new(nu: 32, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu64mum1 = new(nu: 64, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu128mum1 = new(nu: 128, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu129mum1 = new(nu: 129, mu: -1);

        readonly NoncentralStudentTDistribution dist_nu1mu2 = new(nu: 1, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu2mu2 = new(nu: 2, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu3mu2 = new(nu: 3, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu4mu2 = new(nu: 4, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu5mu2 = new(nu: 5, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu8mu2 = new(nu: 8, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu16mu2 = new(nu: 16, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu32mu2 = new(nu: 32, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu64mu2 = new(nu: 64, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu128mu2 = new(nu: 128, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu129mu2 = new(nu: 129, mu: 2);

        readonly NoncentralStudentTDistribution dist_nu1mu4 = new(nu: 1, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu2mu4 = new(nu: 2, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu3mu4 = new(nu: 3, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu4mu4 = new(nu: 4, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu5mu4 = new(nu: 5, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu8mu4 = new(nu: 8, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu16mu4 = new(nu: 16, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu32mu4 = new(nu: 32, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu64mu4 = new(nu: 64, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu128mu4 = new(nu: 128, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu129mu4 = new(nu: 129, mu: 4);

        readonly NoncentralStudentTDistribution dist_nu1mu16 = new(nu: 1, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu2mu16 = new(nu: 2, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu3mu16 = new(nu: 3, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu4mu16 = new(nu: 4, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu5mu16 = new(nu: 5, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu8mu16 = new(nu: 8, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu16mu16 = new(nu: 16, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu32mu16 = new(nu: 32, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu64mu16 = new(nu: 64, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu128mu16 = new(nu: 128, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu129mu16 = new(nu: 129, mu: 16);


        NoncentralStudentTDistribution[] Dists => [
            dist_nu1mu0,
            dist_nu2mu0,
            dist_nu3mu0,
            dist_nu4mu0,
            dist_nu5mu0,
            dist_nu8mu0,
            dist_nu16mu0,
            dist_nu32mu0,
            dist_nu64mu0,
            dist_nu128mu0,
            dist_nu129mu0,

            dist_nu1mu1,
            dist_nu2mu1,
            dist_nu3mu1,
            dist_nu4mu1,
            dist_nu5mu1,
            dist_nu8mu1,
            dist_nu16mu1,
            dist_nu32mu1,
            dist_nu64mu1,
            dist_nu128mu1,
            dist_nu129mu1,

            dist_nu1mum1,
            dist_nu2mum1,
            dist_nu3mum1,
            dist_nu4mum1,
            dist_nu5mum1,
            dist_nu8mum1,
            dist_nu16mum1,
            dist_nu32mum1,
            dist_nu64mum1,
            dist_nu128mum1,
            dist_nu129mum1,

            dist_nu1mu2,
            dist_nu2mu2,
            dist_nu3mu2,
            dist_nu4mu2,
            dist_nu5mu2,
            dist_nu8mu2,
            dist_nu16mu2,
            dist_nu32mu2,
            dist_nu64mu2,
            dist_nu128mu2,
            dist_nu129mu2,

            dist_nu1mu4,
            dist_nu2mu4,
            dist_nu3mu4,
            dist_nu4mu4,
            dist_nu5mu4,
            dist_nu8mu4,
            dist_nu16mu4,
            dist_nu32mu4,
            dist_nu64mu4,
            dist_nu128mu4,
            dist_nu129mu4,

            dist_nu1mu16,
            dist_nu2mu16,
            dist_nu3mu16,
            dist_nu4mu16,
            dist_nu5mu16,
            dist_nu8mu16,
            dist_nu16mu16,
            dist_nu32mu16,
            dist_nu64mu16,
            dist_nu128mu16,
            dist_nu129mu16,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Nu={dist.Nu}");
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
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Console.WriteLine(dist.Mode);

                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode - 1e-25), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode + 1e-25), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Variance)) {
                    continue;
                }

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Kurtosis)) {
                    continue;
                }

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
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
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-24);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-24);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (NoncentralStudentTDistribution dist in Dists) {
                // ignore
                if (dist.Nu >= 32) {
                    continue;
                }

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 40000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.08, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
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
            ddouble[] expected_nu1mu1 = [
                1.029884668413974134e-03,
                1.096791850058143405e-03,
                1.170414928904619503e-03,
                1.251680429571833423e-03,
                1.341679707132120355e-03,
                1.441705276414957052e-03,
                1.553296793463055910e-03,
                1.678299719713448193e-03,
                1.818940807247692096e-03,
                1.977926117154365278e-03,
                2.158569546050187644e-03,
                2.364963132613633057e-03,
                2.602205286002616351e-03,
                2.876710379551117736e-03,
                3.196634275434631113e-03,
                3.572467579000201260e-03,
                4.017875613069993822e-03,
                4.550907851305244049e-03,
                5.195771415030056399e-03,
                5.985483922197922836e-03,
                6.965928400400097330e-03,
                8.202198138483276016e-03,
                9.788777856531904673e-03,
                1.186632215119214941e-02,
                1.465007271803889753e-02,
                1.847926422316873568e-02,
                2.390482638885804006e-02,
                3.184595793883755743e-02,
                4.386030838313985381e-02,
                6.254472411089398443e-02,
                9.181152222380588313e-02,
                1.357256910039695563e-01,
                1.930647052601078273e-01,
                2.494917694932792529e-01,
                2.831609877714570955e-01,
                2.850278366317128276e-01,
                2.635559531170010295e-01,
                2.318845903069630121e-01,
                1.990474065690953087e-01,
                1.692567234400259857e-01,
                1.437974534152678396e-01,
                1.226120165438966403e-01,
                1.051609593824428812e-01,
                9.080820267647864152e-02,
                7.896826530620826290e-02,
                6.914888428225915429e-02,
                6.095231165203827961e-02,
                5.406277735069212748e-02,
                4.823121354230080343e-02,
                4.326125914527433752e-02,
                3.899766085109634267e-02,
                3.531708208481956768e-02,
                3.212097322886890705e-02,
                2.933009258408563955e-02,
                2.688030795991225297e-02,
                2.471937663673680688e-02,
                2.280446787942349077e-02,
                2.110024814800645979e-02,
                1.957739333748615904e-02,
                1.821142616126434147e-02,
                1.698180220217723399e-02,
                1.587118711811738581e-02,
                1.486488160201708542e-02,
                1.395036120227603413e-02,
                1.311690594885521997e-02,
            ];
            ddouble[] expected_nu2mum1 = [
                6.950299731914081486e-03,
                7.606378680696425885e-03,
                8.346270910998968354e-03,
                9.183535475885641272e-03,
                1.013439043360631932e-02,
                1.121832461597420741e-02,
                1.245886998717291375e-02,
                1.388458135405656359e-02,
                1.553028476170790172e-02,
                1.743867504923659814e-02,
                1.966236796672374970e-02,
                2.226654420920740410e-02,
                2.533236257077499734e-02,
                2.896136665818788841e-02,
                3.328116015911077724e-02,
                3.845266717755893571e-02,
                4.467929612254971372e-02,
                5.221821923428458800e-02,
                6.139361749175427374e-02,
                7.261081743675436295e-02,
                8.636815409799364351e-02,
                1.032589733362007239e-01,
                1.239473866188853107e-01,
                1.490850921194628309e-01,
                1.791097459405861791e-01,
                2.138313462481669092e-01,
                2.517021367201676130e-01,
                2.887804414748646464e-01,
                3.178517653196498705e-01,
                3.290764594521785602e-01,
                3.140036633981047043e-01,
                2.722355067279770835e-01,
                2.144409712401767198e-01,
                1.561522156880151835e-01,
                1.081727240454869315e-01,
                7.339914265157548856e-02,
                4.987227471819236091e-02,
                3.438838401021548918e-02,
                2.421844073332200115e-02,
                1.745724475393608047e-02,
                1.287620800135569170e-02,
                9.704695148484886874e-03,
                7.460481959011994527e-03,
                5.838712630538264205e-03,
                4.643490431323359091e-03,
                3.746532114322568673e-03,
                3.062167655886517161e-03,
                2.532069117159483577e-03,
                2.115782224187674147e-03,
                1.784754715856723666e-03,
                1.518502523920622080e-03,
                1.302104504594094048e-03,
                1.124536958997865042e-03,
                9.775480797765156472e-04,
                8.548853266001425451e-04,
                7.517571990380389046e-04,
                6.644530649571906301e-04,
                5.900711121434789198e-04,
                5.263212790508135482e-04,
                4.713808529792548586e-04,
                4.237875143328825585e-04,
                3.823593106112232413e-04,
                3.461342069499561986e-04,
                3.143240131797816962e-04,
                2.862789704052345202e-04,
            ];
            ddouble[] expected_nu4mu2 = [
                1.464098223143212837e-06,
                1.708507273104467589e-06,
                2.003258048036320157e-06,
                2.360795574090711884e-06,
                2.797183852291393688e-06,
                3.333325031285999171e-06,
                3.996648720750331967e-06,
                4.823476971851309684e-06,
                5.862371492414143535e-06,
                7.178926187806678830e-06,
                8.862713948032418555e-06,
                1.103748833084679452e-05,
                1.387637461955258161e-05,
                1.762482691116740770e-05,
                2.263587070250574072e-05,
                2.942511701369128171e-05,
                3.875817803430900600e-05,
                5.179220677816893654e-05,
                7.030967627610245404e-05,
                9.711266865126270997e-05,
                1.367025315612311958e-04,
                1.964779975375205727e-04,
                2.888954726943852662e-04,
                4.354504808875557824e-04,
                6.741632088101218384e-04,
                1.073875594470572771e-03,
                1.761777503561851466e-03,
                2.976102450041582994e-03,
                5.164071161326570271e-03,
                9.151226387250372546e-03,
                1.639131435741447707e-02,
                2.921288759892742490e-02,
                5.075073121372972867e-02,
                8.400084840216752446e-02,
                1.297088034170930382e-01,
                1.840985606707040423e-01,
                2.386935044074002810e-01,
                2.833518736036222396e-01,
                3.105903616371419873e-01,
                3.180849262019659895e-01,
                3.082302682110296677e-01,
                2.859840959118334358e-01,
                2.567019041509658317e-01,
                2.248274648483234162e-01,
                1.934528344581032955e-01,
                1.644120210859889297e-01,
                1.385858320785620657e-01,
                1.162246671264756531e-01,
                9.720865340428283541e-02,
                8.122929998993715095e-02,
                6.790424221775287483e-02,
                5.684351104083793654e-02,
                4.768396826871078992e-02,
                4.010432367417442884e-02,
                3.382909687797335668e-02,
                2.862678339662317814e-02,
                2.430535440822863616e-02,
                2.070685380107149809e-02,
                1.770202219996460866e-02,
                1.518538983886607088e-02,
                1.307100706738975257e-02,
                1.128883665960516564e-02,
                9.781761777549229792e-03,
                8.503134560077370632e-03,
                7.414783465188060259e-03,
            ];
            ddouble[] expected_nu5mu4 = [
                8.102311967558731157e-11,
                9.767061118399378480e-11,
                1.184336829407102899e-10,
                1.445125616877202091e-10,
                1.775143055262130791e-10,
                2.196111890393305120e-10,
                2.737692981765429581e-10,
                3.440828497080931767e-10,
                4.362682987038793346e-10,
                5.584070653250612942e-10,
                7.220757803905575360e-10,
                9.440970382485239705e-10,
                1.249295074099922544e-09,
                1.674916258307741215e-09,
                2.277860919830345043e-09,
                3.146772128659576625e-09,
                4.422716574568579350e-09,
                6.335471527504405433e-09,
                9.268955887636456619e-09,
                1.388279269279282105e-08,
                2.134560372194387322e-08,
                3.379887825163621504e-08,
                5.531555650767860154e-08,
                9.396568530452615008e-08,
                1.664738946017703734e-07,
                3.092425108733207019e-07,
                6.058091264250214947e-07,
                1.258926554047312194e-06,
                2.789775529952722348e-06,
                6.614680686205193721e-06,
                1.677398974043871509e-05,
                4.518872528882844831e-05,
                1.273438577372274937e-04,
                3.660508058119314842e-04,
                1.038384631650249301e-03,
                2.804240917962452957e-03,
                6.980290614366382526e-03,
                1.564494773068389119e-02,
                3.119037111912031229e-02,
                5.521142225374906648e-02,
                8.729512361759993178e-02,
                1.246081248020537990e-01,
                1.626554486077373685e-01,
                1.967276790854739499e-01,
                2.232083568141562113e-01,
                2.402186850569097698e-01,
                2.475755240112292310e-01,
                2.463353809868659150e-01,
                2.382285327193892632e-01,
                2.251792791510779379e-01,
                2.089859678967699308e-01,
                1.911544241138093680e-01,
                1.728445585850971744e-01,
                1.548858489323391707e-01,
                1.378265436601558147e-01,
                1.219934363796226218e-01,
                1.075491995595190542e-01,
                9.454130375409579290e-02,
                8.294079301216454592e-02,
                7.267141685442825294e-02,
                6.363057790748816167e-02,
                5.570379707893942378e-02,
                4.877428637079237600e-02,
                4.272896555268299729e-02,
                3.746197620093304204e-02,
            ];

            foreach ((NoncentralStudentTDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1mu1, expected_nu1mu1),
                (dist_nu2mum1, expected_nu2mum1),
                (dist_nu4mu2, expected_nu4mu2),
                (dist_nu5mu4, expected_nu5mu4),
            }) {
                for ((ddouble x, int i) = (-8, 0); i < expecteds.Length; x += 0.25, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-5, $"{dist} pdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_nu1mu1 = [
                8.285944275278844806e-03,
                8.551647631161951590e-03,
                8.834899481073523256e-03,
                9.137491506887970155e-03,
                9.461466953105700239e-03,
                9.809166200423050791e-03,
                1.018328257422739316e-02,
                1.058693116057262892e-02,
                1.102373429035697489e-02,
                1.149792857205097357e-02,
                1.201450004683440875e-02,
                1.257935641963720333e-02,
                1.319954870611316639e-02,
                1.388355952117703573e-02,
                1.464168238681679696e-02,
                1.548652707272094725e-02,
                1.643370206933785632e-02,
                1.750275006865356989e-02,
                1.871845124069008992e-02,
                2.011267150601997455e-02,
                2.172703537780786703e-02,
                2.361687490475965906e-02,
                2.585720233632439469e-02,
                2.855197691477801697e-02,
                3.184888149342872399e-02,
                3.596356714687924672e-02,
                4.122055920482198960e-02,
                4.812388965619218695e-02,
                5.748009172724949378e-02,
                7.060587579555738236e-02,
                8.963240079750928180e-02,
                1.177449979184550521e-01,
                1.586552539314570187e-01,
                2.142304692010732481e-01,
                2.814609687619795686e-01,
                3.531049214080270726e-01,
                4.220200304964913940e-01,
                4.840473345294080687e-01,
                5.378791827857102303e-01,
                5.838328213835481950e-01,
                6.228719645609039901e-01,
                6.560888868716512334e-01,
                6.844893666933324905e-01,
                7.089273865558834986e-01,
                7.301025986999960526e-01,
                7.485796325277378482e-01,
                7.648120465631097487e-01,
                7.791645234759784211e-01,
                7.919314538726720532e-01,
                8.033518176480952677e-01,
                8.136208584739683669e-01,
                8.228991647194379144e-01,
                8.313197187904681762e-01,
                8.389933786958342354e-01,
                8.460131569626876402e-01,
                8.524575778193083098e-01,
                8.583933264714198597e-01,
                8.638773525873504333e-01,
                8.689585509030792743e-01,
                8.736791123600414766e-01,
                8.780756170466763599e-01,
                8.821799235823520879e-01,
                8.860198970549559849e-01,
                8.896200081509800794e-01,
                8.930018289208374593e-01,
            ];
            ddouble[] expected_nu2mum1 = [
                2.891477558073873161e-02,
                3.073274036959947400e-02,
                3.272494335625369954e-02,
                3.491398283757483900e-02,
                3.732616792357863061e-02,
                3.999225682263339576e-02,
                4.294836729515741358e-02,
                4.623710500526430345e-02,
                4.990896892812979324e-02,
                5.402411058961784779e-02,
                5.865454701533093923e-02,
                6.388695745271558146e-02,
                6.982623308244026439e-02,
                7.659999900911226844e-02,
                8.436439026531357410e-02,
                9.331143791061222159e-02,
                1.036785020206603958e-01,
                1.157602580664462144e-01,
                1.299237586081461171e-01,
                1.466269557443729554e-01,
                1.664405735712990164e-01,
                1.900719416351471791e-01,
                2.183865236293990320e-01,
                2.524169350927559785e-01,
                2.933379503876122341e-01,
                3.423668700663962117e-01,
                4.005227070089473207e-01,
                4.681614424873635216e-01,
                5.442541911098444096e-01,
                6.255975099691507557e-01,
                7.065774638212876102e-01,
                7.803401011833410728e-01,
                8.413447460685430368e-01,
                8.875354666371699786e-01,
                9.203042601355727603e-01,
                9.427357072353101941e-01,
                9.579437363233822555e-01,
                9.683403269956802895e-01,
                9.755783337474334127e-01,
                9.807317211351593844e-01,
                9.844873847471865558e-01,
                9.872865044684904223e-01,
                9.894165285397510390e-01,
                9.910683217742932705e-01,
                9.923712478252055824e-01,
                9.934148070232081684e-01,
                9.942621559217615923e-01,
                9.949586978003077675e-01,
                9.955376423977139178e-01,
                9.960236718676896972e-01,
                9.964354025092152511e-01,
                9.967870669856553345e-01,
                9.970896830430685531e-01,
                9.973518782046488784e-01,
                9.975804802483200229e-01,
                9.977809457818926919e-01,
                9.979576752877422630e-01,
                9.981142474770784423e-01,
                9.982535955645112091e-01,
                9.983781412385092091e-01,
                9.984898974730711529e-01,
                9.985905481481349710e-01,
                9.986815102382238241e-01,
                9.987639827766126022e-01,
                9.988389856989791094e-01,
            ];
            ddouble[] expected_nu4mu2 = [
                2.994963404896093875e-06,
                3.390600056064820445e-06,
                3.853403698686572610e-06,
                4.397450302227477650e-06,
                5.040357945449152211e-06,
                5.804335816953225431e-06,
                6.717593287802726097e-06,
                7.816250601531829836e-06,
                9.146955834507921115e-06,
                1.077050754370086766e-05,
                1.276692668147026699e-05,
                1.524264370709885533e-05,
                1.834081478803688964e-05,
                2.225633395356490032e-05,
                2.725800119580392908e-05,
                3.372177395660454560e-05,
                4.218148430627354645e-05,
                5.340746132574680782e-05,
                6.853166542762460789e-05,
                8.924819353622014118e-05,
                1.181458173387262695e-04,
                1.592681040064256037e-04,
                2.190809552418231110e-04,
                3.081835025758265095e-04,
                4.444044910571695595e-04,
                6.585275543973916068e-04,
                1.005064133808488552e-03,
                1.582728379377389138e-03,
                2.573232017473635336e-03,
                4.312961955230142086e-03,
                7.417406369394169419e-03,
                1.297120788766627231e-02,
                2.275013194817921205e-02,
                3.932870386602460577e-02,
                6.580202153431058409e-02,
                1.049209129184353051e-01,
                1.578735272951267143e-01,
                2.234314909762178392e-01,
                2.980823617018307781e-01,
                3.770676561313910380e-01,
                4.556719904154332701e-01,
                5.301493888185710990e-01,
                5.980805389511567594e-01,
                6.582882483551855124e-01,
                7.105394788654576521e-01,
                7.552117981101059696e-01,
                7.930151019430543569e-01,
                8.247944652750345274e-01,
                8.514066974925347253e-01,
                8.736520362965496389e-01,
                8.922425495953661967e-01,
                9.077927705194003272e-01,
                9.208226177098979282e-01,
                9.317662837446448743e-01,
                9.409833580301443234e-01,
                9.487701316844254773e-01,
                9.553700621270214688e-01,
                9.609829733479805514e-01,
                9.657728955543761762e-01,
                9.698746136933179551e-01,
                9.733990672540551170e-01,
                9.764377658587511943e-01,
                9.790663812625814977e-01,
                9.813476606152297110e-01,
                9.833337859040345297e-01,
            ];
            ddouble[] expected_nu5mu4 = [
                1.317274437534907342e-10,
                1.539888368917298872e-10,
                1.809049263837888536e-10,
                2.136474936128519901e-10,
                2.537365122028485682e-10,
                3.031604526755492999e-10,
                3.645442617139933660e-10,
                4.413869191037125763e-10,
                5.384018912576699753e-10,
                6.620078931368548760e-10,
                8.210519113837876568e-10,
                1.027886545010670908e-09,
                1.299998374737410643e-09,
                1.662516399633673195e-09,
                2.152147882317416718e-09,
                2.823470311699779135e-09,
                3.759215669407108076e-09,
                5.087484877819292706e-09,
                7.011280889201998434e-09,
                9.860648186865716980e-09,
                1.418753758919422125e-08,
                2.094428519505697750e-08,
                3.183191665812046267e-08,
                5.000765068707835326e-08,
                8.158628961307812233e-08,
                1.389812345209406885e-07,
                2.487281163254806415e-07,
                4.708173281242622465e-07,
                9.491636653313761940e-07,
                2.050556992427843284e-06,
                4.765414627461161956e-06,
                1.190166392706082377e-05,
                3.167124183311996498e-05,
                8.817206753183903008e-05,
                2.498017542782856259e-04,
                6.966231979400362706e-04,
                1.850459931492432495e-03,
                4.559685472816781340e-03,
                1.024848915241142328e-02,
                2.086704771453860852e-02,
                3.853325150740972288e-02,
                6.495524596806265261e-02,
                1.008990359268638892e-01,
                1.459484282122374688e-01,
                1.986251769412296642e-01,
                2.567585636285069084e-01,
                3.179265889030191694e-01,
                3.798280740204884043e-01,
                4.405215718691692928e-01,
                4.985308754105037821e-01,
                5.528500566802196925e-01,
                6.028883230771788782e-01,
                6.483883369678906661e-01,
                6.893406495247300647e-01,
                7.259068257922292045e-01,
                7.583565574225711314e-01,
                7.870195997316336012e-01,
                8.122510613500127929e-01,
                8.344076623553093874e-01,
                8.538324567877080273e-01,
                8.708457802841574891e-01,
                8.857405891703987733e-01,
                8.987807735532988840e-01,
                9.102013942322757512e-01,
                9.202100918423642417e-01,
            ];

            foreach ((NoncentralStudentTDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1mu1, expected_nu1mu1),
                (dist_nu2mum1, expected_nu2mum1),
                (dist_nu4mu2, expected_nu4mu2),
                (dist_nu5mu4, expected_nu5mu4),
            }) {
                for ((ddouble x, int i) = (-8, 0); i < expecteds.Length; x += 0.25, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-5, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void PDFIntegrationTest() {
            NoncentralStudentTDistribution dist_p = new(5, 7), dist_n = new(5, -7);

            ddouble expected_m2 = "1.024149272776004956615333514865408381464e-15";
            ddouble expected_m10 = "9.352968724374164571816907170797102136878e-20";
            ddouble expected_p2 = "0.00008255480759792435894808309620397674763859";
            ddouble expected_p10 = "0.07265965604895470390419337567812317002502";

            Assert.IsTrue(ddouble.Abs(expected_m2 - dist_p.PDF(-2)) / expected_m2 < 1e-30);
            Assert.AreEqual(dist_p.PDF(-2), dist_n.PDF(2));

            Assert.IsTrue(ddouble.Abs(expected_m10 - dist_p.PDF(-10)) / expected_m10 < 1e-30);
            Assert.AreEqual(dist_p.PDF(-10), dist_n.PDF(10));

            Assert.IsTrue(ddouble.Abs(expected_p2 - dist_p.PDF(2)) / expected_p2 < 1e-30);
            Assert.AreEqual(dist_p.PDF(2), dist_n.PDF(-2));

            Assert.IsTrue(ddouble.Abs(expected_p10 - dist_p.PDF(10)) / expected_p10 < 1e-30);
            Assert.AreEqual(dist_p.PDF(10), dist_n.PDF(-10));
        }

        [TestMethod()]
        public void PDFLargeMuTest() {
            NoncentralStudentTDistribution dist_1_32 = new(1, 32), dist_1_64 = new(1, 64);
            NoncentralStudentTDistribution dist_2_32 = new(2, 32), dist_2_64 = new(2, 64);

            ddouble expected_1_32_24 = "0.01820401012171220383162340026353465206971";
            ddouble expected_1_32_32 = "0.01510841237253426224138495569370237551505";
            ddouble expected_1_32_48 = "0.008868636760298447264235758063105276307518";
            
            ddouble expected_1_64_48 = "0.009109254125516603938567143674265819735269";
            ddouble expected_1_64_64 = "0.007559739389084211756024140239899313681313";
            ddouble expected_1_64_96 = "0.004436162775362552325578775472944419618799";

            ddouble expected_2_32_24 = "0.02500061513785828548980819708998884927665";
            ddouble expected_2_32_32 = "0.02294765644513584484277875185332015188426";
            ddouble expected_2_32_48 = "0.01186413508809938042012553133888608617727";
            
            ddouble expected_2_64_48 = "0.01251471075031821632571734097288875741963";
            ddouble expected_2_64_64 = "0.01149062222422901172058935997043808705977";
            ddouble expected_2_64_96 = "0.005935657022923590527487703553647298323921";

            Assert.IsTrue(ddouble.Abs(expected_1_32_24 - dist_1_32.PDF(24)) / expected_1_32_24 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_32_32 - dist_1_32.PDF(32)) / expected_1_32_32 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_32_48 - dist_1_32.PDF(48)) / expected_1_32_48 < 2e-29);

            Assert.IsTrue(ddouble.Abs(expected_1_64_48 - dist_1_64.PDF(48)) / expected_1_64_48 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_64_64 - dist_1_64.PDF(64)) / expected_1_64_64 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_1_64_96 - dist_1_64.PDF(96)) / expected_1_64_96 < 2e-29);

            Assert.IsTrue(ddouble.Abs(expected_2_32_24 - dist_2_32.PDF(24)) / expected_2_32_24 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_32_32 - dist_2_32.PDF(32)) / expected_2_32_32 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_32_48 - dist_2_32.PDF(48)) / expected_2_32_48 < 2e-29);

            Assert.IsTrue(ddouble.Abs(expected_2_64_48 - dist_2_64.PDF(48)) / expected_2_64_48 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_64_64 - dist_2_64.PDF(64)) / expected_2_64_64 < 2e-29);
            Assert.IsTrue(ddouble.Abs(expected_2_64_96 - dist_2_64.PDF(96)) / expected_2_64_96 < 2e-29);
        }

        [TestMethod()]
        public void CDFLargeMuTest() {
            NoncentralStudentTDistribution dist_1_32 = new(1, 32), dist_2_32 = new(2, 32);

            ddouble expected_1_32_4 = "8.416933929597894762386824205013044551429e-15";    
            ddouble expected_2_32_4 = "1.8530759205047966738852532488567287089941e-25";
            
            Assert.IsTrue(ddouble.Abs(expected_1_32_4 - dist_1_32.CDF(4)) / expected_1_32_4 < 2e-22);
            Assert.IsTrue(ddouble.Abs(expected_2_32_4 - dist_2_32.CDF(4)) / expected_2_32_4 < 2e-22);
        }
    }
}