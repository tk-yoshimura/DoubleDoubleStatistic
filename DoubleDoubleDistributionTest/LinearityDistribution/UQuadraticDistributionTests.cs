using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.LinearityDistribution {
    [TestClass()]
    public class UQuadraticDistributionTests {
        readonly UQuadraticDistribution dist_a0b1 = new(a: 0, b: 1);
        readonly UQuadraticDistribution dist_a1b2 = new(a: 1, b: 2);
        readonly UQuadraticDistribution dist_a2b4 = new(a: 2, b: 4);

        UQuadraticDistribution[] Dists => [
            dist_a0b1,
            dist_a1b2,
            dist_a2b4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (UQuadraticDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Min={dist.A}");
                Console.WriteLine($"Max={dist.B}");
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
            foreach (UQuadraticDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (UQuadraticDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (UQuadraticDistribution dist in Dists) {
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
            foreach (UQuadraticDistribution dist in Dists) {
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
            foreach (UQuadraticDistribution dist in Dists) {
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
            ddouble[] expected_dist_a0b1 = [
                3.0,
                2.296875,
                1.6875,
                1.171875,
                0.75,
                0.421875,
                0.1875,
                0.046875,
                0.0,
                0.046875,
                0.1875,
                0.421875,
                0.75,
                1.171875,
                1.6875,
                2.296875,
                3.0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ];
            ddouble[] expected_dist_a1b2 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                3.0,
                2.296875,
                1.6875,
                1.171875,
                0.75,
                0.421875,
                0.1875,
                0.046875,
                0.0,
                0.046875,
                0.1875,
                0.421875,
                0.75,
                1.171875,
                1.6875,
                2.296875,
                3.0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ];
            ddouble[] expected_dist_a2b4 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1.5,
                1.318359375,
                1.1484375,
                0.990234375,
                0.84375,
                0.708984375,
                0.5859375,
                0.474609375,
                0.375,
                0.287109375,
                0.2109375,
                0.146484375,
                0.09375,
                0.052734375,
                0.0234375,
                0.005859375,
                0.0,
                0.005859375,
                0.0234375,
                0.052734375,
                0.09375,
                0.146484375,
                0.2109375,
                0.287109375,
                0.375,
                0.474609375,
                0.5859375,
                0.708984375,
                0.84375,
                0.990234375,
                1.1484375,
                1.318359375,
                1.5,
            ];

            foreach ((UQuadraticDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1, expected_dist_a0b1),
                (dist_a1b2, expected_dist_a1b2),
                (dist_a2b4, expected_dist_a2b4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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
            ddouble[] expected_dist_a0b1 = [
                0.0,
                0.1650390625,
                0.2890625,
                0.3779296875,
                0.4374999999999999,
                0.4736328125,
                0.4921875,
                0.4990234374999999,
                0.5,
                0.5009765625,
                0.5078125,
                0.5263671875,
                0.5625,
                0.6220703125,
                0.7109375,
                0.8349609375,
                1.0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ];
            ddouble[] expected_dist_a1b2 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.1650390624999991,
                0.2890625000000009,
                0.3779296875,
                0.4374999999999991,
                0.4736328125000009,
                0.4921875,
                0.4990234374999991,
                0.5000000000000009,
                0.5009765625,
                0.5078124999999991,
                0.5263671875000009,
                0.5625,
                0.6220703124999991,
                0.7109375000000009,
                0.8349609375,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ];
            ddouble[] expected_dist_a2b4 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.08801269531250089,
                0.1650390624999991,
                0.2318115234375,
                0.2890625000000009,
                0.3375244140624991,
                0.3779296875,
                0.4110107421875009,
                0.4374999999999991,
                0.4581298828125,
                0.4736328125000009,
                0.4847412109374991,
                0.4921875,
                0.4967041015625009,
                0.4990234374999991,
                0.4998779296875,
                0.5000000000000009,
                0.5001220703124991,
                0.5009765625,
                0.5032958984375009,
                0.5078124999999991,
                0.5152587890625,
                0.5263671875000009,
                0.5418701171874991,
                0.5625,
                0.5889892578125009,
                0.6220703124999991,
                0.6624755859375,
                0.7109375000000009,
                0.7681884765624991,
                0.8349609375,
                0.9119873046875009,
                1,
            ];

            foreach ((UQuadraticDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1, expected_dist_a0b1),
                (dist_a1b2, expected_dist_a1b2),
                (dist_a2b4, expected_dist_a2b4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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