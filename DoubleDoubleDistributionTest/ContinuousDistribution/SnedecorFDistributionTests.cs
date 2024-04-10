using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ContinuousDistribution {
    [TestClass()]
    public class SnedecorFDistributionTests {
        readonly SnedecorFDistribution dist_d1_1_d2_1 = new(d1: 1, d2: 1);
        readonly SnedecorFDistribution dist_d1_2_d2_1 = new(d1: 2, d2: 1);
        readonly SnedecorFDistribution dist_d1_1_d2_2 = new(d1: 1, d2: 2);
        readonly SnedecorFDistribution dist_d1_2_d2_2 = new(d1: 2, d2: 2);
        readonly SnedecorFDistribution dist_d1_3_d2_4 = new(d1: 3, d2: 4);
        readonly SnedecorFDistribution dist_d1_4_d2_2 = new(d1: 4, d2: 2);
        readonly SnedecorFDistribution dist_d1_2_d2_4 = new(d1: 2, d2: 4);
        readonly SnedecorFDistribution dist_d1_4_d2_4 = new(d1: 4, d2: 4);
        readonly SnedecorFDistribution dist_d1_6_d2_8 = new(d1: 6, d2: 8);
        readonly SnedecorFDistribution dist_d1_8_d2_6 = new(d1: 8, d2: 6);
        readonly SnedecorFDistribution dist_d1_8_d2_8 = new(d1: 8, d2: 8);
        readonly SnedecorFDistribution dist_d1_10_d2_10 = new(d1: 10, d2: 10);

        SnedecorFDistribution[] Dists => [
            dist_d1_1_d2_1,
            dist_d1_2_d2_1,
            dist_d1_1_d2_2,
            dist_d1_2_d2_2,
            dist_d1_3_d2_4,
            dist_d1_4_d2_2,
            dist_d1_2_d2_4,
            dist_d1_4_d2_4,
            dist_d1_6_d2_8,
            dist_d1_8_d2_6,
            dist_d1_8_d2_8,
            dist_d1_10_d2_10
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"D1={dist.D1}");
                Console.WriteLine($"D2={dist.D2}");
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
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (SnedecorFDistribution dist in Dists) {
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
            foreach (SnedecorFDistribution dist in Dists) {
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
            foreach (SnedecorFDistribution dist in Dists) {
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
            ddouble[] expected_dist_d1_1_d2_1 = [
            ];
            ddouble[] expected_dist_d1_2_d2_1 = [
            ];
            ddouble[] expected_dist_d1_1_d2_2 = [
            ];
            ddouble[] expected_dist_d1_2_d2_2 = [
            ];
            ddouble[] expected_dist_d1_3_d2_4 = [
            ];
            ddouble[] expected_dist_d1_4_d2_2 = [
            ];
            ddouble[] expected_dist_d1_2_d2_4 = [
            ];
            ddouble[] expected_dist_d1_4_d2_4 = [
            ];
            ddouble[] expected_dist_d1_6_d2_8 = [
            ];
            ddouble[] expected_dist_d1_8_d2_6 = [
            ];
            ddouble[] expected_dist_d1_8_d2_8 = [
            ];
            ddouble[] expected_dist_d1_10_d2_10 = [
            ];

            foreach ((SnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_d1_1_d2_1,   expected_dist_d1_1_d2_1),
                 (dist_d1_2_d2_1,   expected_dist_d1_2_d2_1),
                 (dist_d1_1_d2_2,   expected_dist_d1_1_d2_2),
                 (dist_d1_2_d2_2,   expected_dist_d1_2_d2_2),
                 (dist_d1_3_d2_4,   expected_dist_d1_3_d2_4),
                 (dist_d1_4_d2_2,   expected_dist_d1_4_d2_2),
                 (dist_d1_2_d2_4,   expected_dist_d1_2_d2_4),
                 (dist_d1_4_d2_4,   expected_dist_d1_4_d2_4),
                 (dist_d1_6_d2_8,   expected_dist_d1_6_d2_8),
                 (dist_d1_8_d2_6,   expected_dist_d1_8_d2_6),
                 (dist_d1_8_d2_8,   expected_dist_d1_8_d2_8),
                 (dist_d1_10_d2_10, expected_dist_d1_10_d2_10),
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
            ddouble[] expected_dist_d1_1_d2_1 = [
            ];
            ddouble[] expected_dist_d1_2_d2_1 = [
            ];
            ddouble[] expected_dist_d1_1_d2_2 = [
            ];
            ddouble[] expected_dist_d1_2_d2_2 = [
            ];
            ddouble[] expected_dist_d1_3_d2_4 = [
            ];
            ddouble[] expected_dist_d1_4_d2_2 = [
            ];
            ddouble[] expected_dist_d1_2_d2_4 = [
            ];
            ddouble[] expected_dist_d1_4_d2_4 = [
            ];
            ddouble[] expected_dist_d1_6_d2_8 = [
            ];
            ddouble[] expected_dist_d1_8_d2_6 = [
            ];
            ddouble[] expected_dist_d1_8_d2_8 = [
            ];
            ddouble[] expected_dist_d1_10_d2_10 = [
            ];

            foreach ((SnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_d1_1_d2_1,   expected_dist_d1_1_d2_1),
                 (dist_d1_2_d2_1,   expected_dist_d1_2_d2_1),
                 (dist_d1_1_d2_2,   expected_dist_d1_1_d2_2),
                 (dist_d1_2_d2_2,   expected_dist_d1_2_d2_2),
                 (dist_d1_3_d2_4,   expected_dist_d1_3_d2_4),
                 (dist_d1_4_d2_2,   expected_dist_d1_4_d2_2),
                 (dist_d1_2_d2_4,   expected_dist_d1_2_d2_4),
                 (dist_d1_4_d2_4,   expected_dist_d1_4_d2_4),
                 (dist_d1_6_d2_8,   expected_dist_d1_6_d2_8),
                 (dist_d1_8_d2_6,   expected_dist_d1_8_d2_6),
                 (dist_d1_8_d2_8,   expected_dist_d1_8_d2_8),
                 (dist_d1_10_d2_10, expected_dist_d1_10_d2_10),
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