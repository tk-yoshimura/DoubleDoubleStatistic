using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ContinuousDistribution {
    [TestClass()]
    public class IrwinHallDistributionTests {
        readonly IrwinHallDistribution dist_n1 = new(n: 1);
        readonly IrwinHallDistribution dist_n2 = new(n: 2);
        readonly IrwinHallDistribution dist_n3 = new(n: 3);
        readonly IrwinHallDistribution dist_n4 = new(n: 4);
        readonly IrwinHallDistribution dist_n5 = new(n: 5);
        readonly IrwinHallDistribution dist_n8 = new(n: 8);
        readonly IrwinHallDistribution dist_n16 = new(n: 16);
        readonly IrwinHallDistribution dist_n32 = new(n: 32);
        readonly IrwinHallDistribution dist_n64 = new(n: 64);


        IrwinHallDistribution[] Dists => new[]{
            dist_n1,
            dist_n2,
            dist_n3,
            dist_n4,
            dist_n5,
            dist_n8,
            dist_n16,
            dist_n32,
            dist_n64,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"N={dist.N}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
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
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= dist.N; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);

                for (ddouble x = 0; x <= dist.N; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);

                for (ddouble x = 0; x <= dist.N; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-22);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 100; i++) {
                    ddouble p = (ddouble)i / 100;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-20);
                    }

                    Assert.IsTrue(ddouble.IsFinite(x));
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 100; i++) {
                    ddouble p = (ddouble)i / 100;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-20);
                    }

                    Assert.IsTrue(ddouble.IsFinite(x));
                }
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_n1 = [
            ];
            ddouble[] expected_n2 = [
            ];
            ddouble[] expected_n4 = [
            ];

            foreach ((IrwinHallDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n1, expected_n1), (dist_n2, expected_n2), (dist_n4, expected_n4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
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
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_n1 = [
            ];
            ddouble[] expected_n2 = [
            ];
            ddouble[] expected_n4 = [
            ];

            foreach ((IrwinHallDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n1, expected_n1), (dist_n2, expected_n2), (dist_n4, expected_n4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }
    }
}