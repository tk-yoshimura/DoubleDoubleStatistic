using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ContinuousDistribution {
    [TestClass()]
    public class GammaDistributionTests {
        readonly GammaDistribution dist_k1theta1 = new(k: 1, theta: 1);
        readonly GammaDistribution dist_k2theta1 = new(k: 2, theta: 1);
        readonly GammaDistribution dist_k1theta2 = new(k: 1, theta: 2);
        readonly GammaDistribution dist_k2theta2 = new(k: 2, theta: 2);
        readonly GammaDistribution dist_k3theta4 = new(k: 3, theta: 4);

        GammaDistribution[] Dists => new[]{
            dist_k1theta1,
            dist_k2theta1,
            dist_k1theta2,
            dist_k2theta2,
            dist_k3theta4,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.K}");
                Console.WriteLine($"Theta={dist.Theta}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            ddouble[] expected_dist_k1theta1 = [
            ];
            ddouble[] expected_dist_k2theta1 = [
            ];
            ddouble[] expected_dist_k1theta2 = [
            ];
            ddouble[] expected_dist_k2theta2 = [
            ];
            ddouble[] expected_dist_k3theta4 = [
            ];

            foreach ((GammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1theta1, expected_dist_k1theta1), (dist_k2theta1, expected_dist_k2theta1),
                (dist_k1theta2, expected_dist_k1theta2), (dist_k2theta2, expected_dist_k2theta2),
                (dist_k3theta4, expected_dist_k3theta4),
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
            ddouble[] expected_dist_k1theta1 = [
            ];
            ddouble[] expected_dist_k2theta1 = [
            ];
            ddouble[] expected_dist_k1theta2 = [
            ];
            ddouble[] expected_dist_k2theta2 = [
            ];
            ddouble[] expected_dist_k3theta4 = [
            ];

            foreach ((GammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1theta1, expected_dist_k1theta1), (dist_k2theta1, expected_dist_k2theta1),
                (dist_k1theta2, expected_dist_k1theta2), (dist_k2theta2, expected_dist_k2theta2),
                (dist_k3theta4, expected_dist_k3theta4),
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