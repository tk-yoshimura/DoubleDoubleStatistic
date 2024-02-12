using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ContinuousDistribution {
    [TestClass()]
    public class ParetoDistributionTests {
        readonly ParetoDistribution dist_xm1alpha2 = new(xm: 1, alpha: 2);
        readonly ParetoDistribution dist_xm2alpha3 = new(xm: 2, alpha: 3);
        readonly ParetoDistribution dist_xm2alpha4 = new(xm: 2, alpha: 4);
        readonly ParetoDistribution dist_xm3alpha4 = new(xm: 3, alpha: 4);

        ParetoDistribution[] Dists => new[]{
            dist_xm1alpha2,
            dist_xm2alpha3,
            dist_xm2alpha4,
            dist_xm3alpha4,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Xm={dist.Xm}");
                Console.WriteLine($"Alpha={dist.Alpha}");
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
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
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
            ddouble[] expected_dist_xm1alpha2 = [
            ];
            ddouble[] expected_dist_xm2alpha3 = [
            ];
            ddouble[] expected_dist_xm2alpha4 = [
            ];
            ddouble[] expected_dist_xm3alpha4 = [
            ];

            foreach ((ParetoDistribution dist, ddouble[] expecteds) in new[]{
                (dist_xm1alpha2, expected_dist_xm1alpha2),
                (dist_xm2alpha3, expected_dist_xm2alpha3),
                (dist_xm2alpha4, expected_dist_xm2alpha4),
                (dist_xm3alpha4, expected_dist_xm3alpha4),
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
            ddouble[] expected_dist_xm1alpha2 = [
            ];
            ddouble[] expected_dist_xm2alpha3 = [
            ];
            ddouble[] expected_dist_xm2alpha4 = [
            ];
            ddouble[] expected_dist_xm3alpha4 = [
            ];

            foreach ((ParetoDistribution dist, ddouble[] expecteds) in new[]{
                (dist_xm1alpha2, expected_dist_xm1alpha2),
                (dist_xm2alpha3, expected_dist_xm2alpha3),
                (dist_xm2alpha4, expected_dist_xm2alpha4),
                (dist_xm3alpha4, expected_dist_xm3alpha4),
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