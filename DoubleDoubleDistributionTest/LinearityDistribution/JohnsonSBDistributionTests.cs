
using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.LinearityDistribution {
    [TestClass()]
    public class JohnsonSBDistributionTests {
        readonly JohnsonSBDistribution dist_gamma1delta1mu0sigma1 = new(gamma: 1, delta: 1, mu: 0, sigma: 1);
        readonly JohnsonSBDistribution dist_gamma2delta1mu1sigma1 = new(gamma: 2, delta: 1, mu: 1, sigma: 1);
        readonly JohnsonSBDistribution dist_gamma1delta1mu0sigma2 = new(gamma: 1, delta: 1, mu: 0, sigma: 2);
        readonly JohnsonSBDistribution dist_gamma2delta2mu2sigma2 = new(gamma: 2, delta: 2, mu: 2, sigma: 2);
        readonly JohnsonSBDistribution dist_gamma3delta2mu0sigma4 = new(gamma: 3, delta: 2, mu: 0, sigma: 4);

        JohnsonSBDistribution[] Dists => [
            dist_gamma1delta1mu0sigma1,
            dist_gamma2delta1mu1sigma1,
            dist_gamma1delta1mu0sigma2,
            dist_gamma2delta2mu2sigma2,
            dist_gamma3delta2mu0sigma4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Gamma={dist.Gamma}");
                Console.WriteLine($"Delta={dist.Delta}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Sigma={dist.Sigma}");
                Console.WriteLine($"Median={dist.Median}");

                /* TODO: Implement */
                //Console.WriteLine($"Mean={dist.Mean}");
                //Console.WriteLine($"Mode={dist.Mode}");
                //Console.WriteLine($"Variance={dist.Variance}");
                //Console.WriteLine($"Skewness={dist.Skewness}");
                //Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                //Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            ddouble[] expected_dist_gamma1delta1mu0sigma1 = [];
            ddouble[] expected_dist_gamma2delta1mu1sigma1 = [];
            ddouble[] expected_dist_gamma1delta1mu0sigma2 = [];
            ddouble[] expected_dist_gamma2delta2mu2sigma2 = [];
            ddouble[] expected_dist_gamma3delta2mu0sigma4 = [];

            foreach ((JohnsonSBDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1delta1mu0sigma1, expected_dist_gamma1delta1mu0sigma1),
                (dist_gamma2delta1mu1sigma1, expected_dist_gamma2delta1mu1sigma1),
                (dist_gamma1delta1mu0sigma2, expected_dist_gamma1delta1mu0sigma2),
                (dist_gamma2delta2mu2sigma2, expected_dist_gamma2delta2mu2sigma2),
                (dist_gamma3delta2mu0sigma4, expected_dist_gamma3delta2mu0sigma4),
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
            ddouble[] expected_dist_gamma1delta1mu0sigma1 = [];
            ddouble[] expected_dist_gamma2delta1mu1sigma1 = [];
            ddouble[] expected_dist_gamma1delta1mu0sigma2 = [];
            ddouble[] expected_dist_gamma2delta2mu2sigma2 = [];
            ddouble[] expected_dist_gamma3delta2mu0sigma4 = [];

            foreach ((JohnsonSBDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1delta1mu0sigma1, expected_dist_gamma1delta1mu0sigma1),
                (dist_gamma2delta1mu1sigma1, expected_dist_gamma2delta1mu1sigma1),
                (dist_gamma1delta1mu0sigma2, expected_dist_gamma1delta1mu0sigma2),
                (dist_gamma2delta2mu2sigma2, expected_dist_gamma2delta2mu2sigma2),
                (dist_gamma3delta2mu0sigma4, expected_dist_gamma3delta2mu0sigma4),
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