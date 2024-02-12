using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ContinuousDistribution {
    [TestClass()]
    public class ReciprocalDistributionTests {
        readonly ReciprocalDistribution dist_a1b2 = new(a: 1, b: 2);
        readonly ReciprocalDistribution dist_a2b3 = new(a: 2, b: 3);
        readonly ReciprocalDistribution dist_a2b4 = new(a: 2, b: 4);
        readonly ReciprocalDistribution dist_a3b4 = new(a: 3, b: 4);

        ReciprocalDistribution[] Dists => new[]{
            dist_a1b2,
            dist_a2b3,
            dist_a2b4,
            dist_a3b4,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (ReciprocalDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"A={dist.A}");
                Console.WriteLine($"B={dist.B}");
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
            foreach (ReciprocalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ReciprocalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ReciprocalDistribution dist in Dists) {
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
            foreach (ReciprocalDistribution dist in Dists) {
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
            foreach (ReciprocalDistribution dist in Dists) {
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
        public void InfoExpectedTest() {
            Assert.IsTrue(ddouble.Abs("1.442695040888963407359924681001892137427" - dist_a1b2.Mean) < 1e-26);
            Assert.IsTrue(ddouble.Abs("8.26735803278373131703054177678467810759e-2" - dist_a1b2.Variance) < 1e-26);
            Assert.IsTrue(ddouble.Abs("2.393417858997472994932677036954518839165e-1" - dist_a1b2.Skewness) < 1e-26);
            Assert.IsTrue(ddouble.Abs("-1.137246458476132118678312005465854454344" - dist_a1b2.Kurtosis) < 1e-26);
            Assert.IsTrue(ddouble.Abs("-1.993933030169167230382309750358118541653e-2" - dist_a1b2.Entropy) < 1e-26);

            Assert.IsTrue(ddouble.Abs("2.885390081777926814719849362003784274853" - dist_a2b4.Mean) < 1e-26);
            Assert.IsTrue(ddouble.Abs("3.306943213113492526812216710713871243047e-1" - dist_a2b4.Variance) < 1e-26);
            Assert.IsTrue(ddouble.Abs("2.393417858997472994932677036954518839265e-1" - dist_a2b4.Skewness) < 1e-26);
            Assert.IsTrue(ddouble.Abs("-1.137246458476132118678312005465854454344" - dist_a2b4.Kurtosis) < 1e-26);
            Assert.IsTrue(ddouble.Abs("6.73207850258253637113409023954595382659e-1" - dist_a2b4.Entropy) < 1e-26);

            Assert.IsTrue(ddouble.Abs("3.476059496782206910376499401645766368748" - dist_a3b4.Mean) < 1e-26);
            Assert.IsTrue(ddouble.Abs("8.321861356795465240379689446126489475011e-2" - dist_a3b4.Variance) < 1e-26);
            Assert.IsTrue(ddouble.Abs("9.960043892185942044437912059971983108507e-2" - dist_a3b4.Skewness) < 1e-26);
            Assert.IsTrue(ddouble.Abs("-1.189134564437105737117615891768726318608" - dist_a3b4.Kurtosis) < 1e-26);
            Assert.IsTrue(ddouble.Abs("-3.445998813238043265926062980199643557012e-3" - dist_a3b4.Entropy) < 1e-26);
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist_a1b2 = [
            ];
            ddouble[] expected_dist_a2b3 = [
            ];
            ddouble[] expected_dist_a2b4 = [
            ];
            ddouble[] expected_dist_a3b4 = [
            ];

            foreach ((ReciprocalDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b2, expected_dist_a1b2),
                (dist_a2b3, expected_dist_a2b3),
                (dist_a2b4, expected_dist_a2b4),
                (dist_a3b4, expected_dist_a3b4),
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
            ddouble[] expected_dist_a1b2 = [
            ];
            ddouble[] expected_dist_a2b3 = [
            ];
            ddouble[] expected_dist_a2b4 = [
            ];
            ddouble[] expected_dist_a3b4 = [
            ];

            foreach ((ReciprocalDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b2, expected_dist_a1b2),
                (dist_a2b3, expected_dist_a2b3),
                (dist_a2b4, expected_dist_a2b4),
                (dist_a3b4, expected_dist_a3b4),
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