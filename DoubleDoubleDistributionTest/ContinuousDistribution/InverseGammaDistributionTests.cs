using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ContinuousDistribution {
    [TestClass()]
    public class InverseGammaDistributionTests {
        readonly InverseGammaDistribution dist_alpha1beta1 = new(alpha: 1, beta: 1);
        readonly InverseGammaDistribution dist_alpha2beta1 = new(alpha: 2, beta: 1);
        readonly InverseGammaDistribution dist_alpha1beta2 = new(alpha: 1, beta: 2);
        readonly InverseGammaDistribution dist_alpha2beta2 = new(alpha: 2, beta: 2);
        readonly InverseGammaDistribution dist_alpha3beta4 = new(alpha: 3, beta: 4);

        InverseGammaDistribution[] Dists => [
            dist_alpha1beta1,
            dist_alpha2beta1,
            dist_alpha1beta2,
            dist_alpha2beta2,
            dist_alpha3beta4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Beta={dist.Beta}");
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
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha1beta1 = [
            ];
            ddouble[] expected_dist_alpha2beta1 = [
            ];
            ddouble[] expected_dist_alpha1beta2 = [
            ];
            ddouble[] expected_dist_alpha2beta2 = [
            ];
            ddouble[] expected_dist_alpha3beta4 = [
            ];

            foreach ((InverseGammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1beta1, expected_dist_alpha1beta1), (dist_alpha2beta1, expected_dist_alpha2beta1),
                (dist_alpha1beta2, expected_dist_alpha1beta2), (dist_alpha2beta2, expected_dist_alpha2beta2),
                (dist_alpha3beta4, expected_dist_alpha3beta4),
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
            ddouble[] expected_dist_alpha1beta1 = [
            ];
            ddouble[] expected_dist_alpha2beta1 = [
            ];
            ddouble[] expected_dist_alpha1beta2 = [
            ];
            ddouble[] expected_dist_alpha2beta2 = [
            ];
            ddouble[] expected_dist_alpha3beta4 = [
            ];

            foreach ((InverseGammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1beta1, expected_dist_alpha1beta1), (dist_alpha2beta1, expected_dist_alpha2beta1),
                (dist_alpha1beta2, expected_dist_alpha1beta2), (dist_alpha2beta2, expected_dist_alpha2beta2),
                (dist_alpha3beta4, expected_dist_alpha3beta4),
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