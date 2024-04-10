using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.LinearityDistribution {
    [TestClass()]
    public class SkewNormalDistributionTests {
        readonly SkewNormalDistribution dist_alpha1mu1sigma1 = new(alpha: 1, mu: 1, sigma: 1);
        readonly SkewNormalDistribution dist_alpha2mu1sigma2 = new(alpha: 2, mu: 1, sigma: 2);
        readonly SkewNormalDistribution dist_alpha1mu2sigma1 = new(alpha: 1, mu: 2, sigma: 1);
        readonly SkewNormalDistribution dist_alpha2mu2sigma2 = new(alpha: 2, mu: 2, sigma: 2);
        readonly SkewNormalDistribution dist_alpha3mu4sigma3 = new(alpha: 3, mu: 4, sigma: 3);

        SkewNormalDistribution[] Dists => [
            dist_alpha1mu1sigma1,
            dist_alpha2mu1sigma2,
            dist_alpha1mu2sigma1,
            dist_alpha2mu2sigma2,
            dist_alpha3mu4sigma3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (SkewNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Sigma={dist.Sigma}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                /* TODO: implement entropy */
                //Console.WriteLine($"Median={dist.Median}");
                //Console.WriteLine($"Mode={dist.Mode}");
                //Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (SkewNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (SkewNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (SkewNormalDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (SkewNormalDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (SkewNormalDistribution dist in Dists) {
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
        public void PDFExpectedTest() {
            ddouble[] expected_dist_alpha1mu1sigma1 = [
            ];
            ddouble[] expected_dist_alpha2mu1sigma2 = [
            ];
            ddouble[] expected_dist_alpha1mu2sigma1 = [
            ];
            ddouble[] expected_dist_alpha2mu2sigma2 = [
            ];
            ddouble[] expected_dist_alpha3mu4sigma3 = [
            ];

            foreach ((SkewNormalDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1mu1sigma1, expected_dist_alpha1mu1sigma1),
                (dist_alpha2mu1sigma2, expected_dist_alpha2mu1sigma2),
                (dist_alpha1mu2sigma1, expected_dist_alpha1mu2sigma1),
                (dist_alpha2mu2sigma2, expected_dist_alpha2mu2sigma2),
                (dist_alpha3mu4sigma3, expected_dist_alpha3mu4sigma3),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 1d / 64, i++) {
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
            ddouble[] expected_dist_alpha1mu1sigma1 = [
            ];
            ddouble[] expected_dist_alpha2mu1sigma2 = [
            ];
            ddouble[] expected_dist_alpha1mu2sigma1 = [
            ];
            ddouble[] expected_dist_alpha2mu2sigma2 = [
            ];
            ddouble[] expected_dist_alpha3mu4sigma3 = [
            ];

            foreach ((SkewNormalDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1mu1sigma1, expected_dist_alpha1mu1sigma1),
                (dist_alpha2mu1sigma2, expected_dist_alpha2mu1sigma2),
                (dist_alpha1mu2sigma1, expected_dist_alpha1mu2sigma1),
                (dist_alpha2mu2sigma2, expected_dist_alpha2mu2sigma2),
                (dist_alpha3mu4sigma3, expected_dist_alpha3mu4sigma3),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 1d / 64, i++) {
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