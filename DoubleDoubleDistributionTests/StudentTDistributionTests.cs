using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistribution.Tests {
    [TestClass()]
    public class StudentTDistributionTests {
        readonly StudentTDistribution dist_nu0p5 = new(nu: 0.5);
        readonly StudentTDistribution dist_nu1 = new(nu: 1);
        readonly StudentTDistribution dist_nu1p5 = new(nu: 1.5);
        readonly StudentTDistribution dist_nu2 = new(nu: 2);
        readonly StudentTDistribution dist_nu2p5 = new(nu: 2.5);
        readonly StudentTDistribution dist_nu3 = new(nu: 3);
        readonly StudentTDistribution dist_nu3p5 = new(nu: 3.5);
        readonly StudentTDistribution dist_nu4 = new(nu: 4);
        readonly StudentTDistribution dist_nu4p5 = new(nu: 4.5);
        readonly StudentTDistribution dist_nu5 = new(nu: 5);

        StudentTDistribution[] Dists => new[]{
            dist_nu0p5,
            dist_nu1,
            dist_nu1p5,
            dist_nu2,
            dist_nu2p5,
            dist_nu3,
            dist_nu3p5,
            dist_nu4,
            dist_nu4p5,
            dist_nu5,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mu={dist.Nu}");
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
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    Console.WriteLine($"pdf({x})={dist.PDF(x)}");
                }
            }
        }

        [TestMethod()]
        public void CDFTest() {
            Assert.Inconclusive();

            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    Console.WriteLine($"cdf({x})={dist.CDF(x)}");
                }
            }
        }

        [TestMethod()]
        public void QuantileTest() {
            Assert.Inconclusive();

            foreach (StudentTDistribution dist in Dists) {
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