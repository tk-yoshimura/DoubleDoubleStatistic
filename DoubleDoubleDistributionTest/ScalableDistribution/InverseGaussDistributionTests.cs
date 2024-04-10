using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ScalableDistribution {
    [TestClass()]
    public class InverseGaussDistributionTests {
        readonly InverseGaussDistribution dist_mu1lambda1 = new(mu: 1, lambda: 1);
        readonly InverseGaussDistribution dist_mu2lambda1 = new(mu: 2, lambda: 1);
        readonly InverseGaussDistribution dist_mu1lambda2 = new(mu: 1, lambda: 2);
        readonly InverseGaussDistribution dist_mu2lambda2 = new(mu: 2, lambda: 2);
        readonly InverseGaussDistribution dist_mu3lambda4 = new(mu: 3, lambda: 4);

        InverseGaussDistribution[] Dists => [
            dist_mu1lambda1,
            dist_mu2lambda1,
            dist_mu1lambda2,
            dist_mu2lambda2,
            dist_mu3lambda4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Lambda={dist.Lambda}");
                Console.WriteLine($"Mean={dist.Mean}");
                /* TODO: Implement */
                //Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                /* TODO: Implement */
                //Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (InverseGaussDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (InverseGaussDistribution dist in Dists) {
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

            foreach (InverseGaussDistribution dist in Dists) {
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
            ddouble[] expected_dist_mu1lambda1 = [
            ];
            ddouble[] expected_dist_mu2lambda1 = [
            ];
            ddouble[] expected_dist_mu1lambda2 = [
            ];
            ddouble[] expected_dist_mu2lambda2 = [
            ];
            ddouble[] expected_dist_mu3lambda4 = [
            ];

            foreach ((InverseGaussDistribution dist, ddouble[] expecteds) in new[]{
                (dist_mu1lambda1, expected_dist_mu1lambda1),
                (dist_mu2lambda1, expected_dist_mu2lambda1),
                (dist_mu1lambda2, expected_dist_mu1lambda2),
                (dist_mu2lambda2, expected_dist_mu2lambda2),
                (dist_mu3lambda4, expected_dist_mu3lambda4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
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
            ddouble[] expected_dist_mu1lambda1 = [
            ];
            ddouble[] expected_dist_mu2lambda1 = [
            ];
            ddouble[] expected_dist_mu1lambda2 = [
            ];
            ddouble[] expected_dist_mu2lambda2 = [
            ];
            ddouble[] expected_dist_mu3lambda4 = [
            ];

            foreach ((InverseGaussDistribution dist, ddouble[] expecteds) in new[]{
                (dist_mu1lambda1, expected_dist_mu1lambda1),
                (dist_mu2lambda1, expected_dist_mu2lambda1),
                (dist_mu1lambda2, expected_dist_mu1lambda2),
                (dist_mu2lambda2, expected_dist_mu2lambda2),
                (dist_mu3lambda4, expected_dist_mu3lambda4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }
    }
}