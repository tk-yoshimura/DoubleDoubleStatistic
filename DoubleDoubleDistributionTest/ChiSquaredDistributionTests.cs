using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistribution.Tests {
    [TestClass()]
    public class ChiSquaredDistributionTests {
        readonly ChiSquaredDistribution dist_nu1 = new(k: 1);
        readonly ChiSquaredDistribution dist_nu2 = new(k: 2);
        readonly ChiSquaredDistribution dist_nu3 = new(k: 3);
        readonly ChiSquaredDistribution dist_nu4 = new(k: 4);
        readonly ChiSquaredDistribution dist_nu5 = new(k: 5);
        readonly ChiSquaredDistribution dist_nu8 = new(k: 8);
        readonly ChiSquaredDistribution dist_nu16 = new(k: 16);
        readonly ChiSquaredDistribution dist_nu32 = new(k: 32);
        readonly ChiSquaredDistribution dist_nu64 = new(k: 64);
        readonly ChiSquaredDistribution dist_nu128 = new(k: 128);


        ChiSquaredDistribution[] Dists => new[]{
            dist_nu1,
            dist_nu2,
            dist_nu3,
            dist_nu4,
            dist_nu5,
            dist_nu8,
            dist_nu16,
            dist_nu32,
            dist_nu64,
            dist_nu128,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (ChiSquaredDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.K}");
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
            foreach (ChiSquaredDistribution dist in Dists) {
                Console.WriteLine(dist);
                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    Console.WriteLine($"pdf({x})={dist.PDF(x)}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    Console.WriteLine($"pdf({x})={dist.PDF(x)}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    Console.WriteLine($"pdf({x})={dist.PDF(x)}");
                }
            }
        }

        [TestMethod()]
        public void CDFTest() {
            foreach (ChiSquaredDistribution dist in Dists) {
                Console.WriteLine(dist);

                Console.WriteLine($"cdf(median)={dist.CDF(dist.Median)}");

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    Console.WriteLine($"cdf({x})={dist.CDF(x)}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    Console.WriteLine($"cdf({x})={dist.CDF(x)}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    Console.WriteLine($"cdf({x})={dist.CDF(x)}");
                }

            }
        }

        [TestMethod()]
        public void QuantileTest() {
            Assert.Inconclusive();

            foreach (ChiSquaredDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={dist.CDF(x)}");
                }
            }
        }
    }
}