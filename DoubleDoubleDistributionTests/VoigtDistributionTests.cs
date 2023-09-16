using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistribution.Tests {
    [TestClass()]
    public class VoigtDistributionTests {
        VoigtDistribution dist1 = new();
        VoigtDistribution dist2 = new(2, 3);

        [TestMethod()]
        public void VoigtDistributionTest() {
            Console.WriteLine(dist1);
            Console.WriteLine($"Domain={dist1.Domain}");
            Console.WriteLine($"Gamma={dist1.Gamma}");
            Console.WriteLine($"Sigma={dist1.Sigma}");
            Console.WriteLine($"Mean={dist1.Mean}");
            Console.WriteLine($"Median={dist1.Median}");
            Console.WriteLine($"Mode={dist1.Mode}");
            Console.WriteLine($"Variance={dist1.Variance}");
            Console.WriteLine($"Skewness={dist1.Skewness}");
            Console.WriteLine($"Kurtosis={dist1.Kurtosis}");
            Console.WriteLine($"Entropy={dist1.Entropy}");

            Console.WriteLine(dist2);
            Console.WriteLine($"Domain={dist2.Domain}");
            Console.WriteLine($"Gamma={dist2.Gamma}");
            Console.WriteLine($"Sigma={dist2.Sigma}");
            Console.WriteLine($"Mean={dist2.Mean}");
            Console.WriteLine($"Median={dist2.Median}");
            Console.WriteLine($"Mode={dist2.Mode}");
            Console.WriteLine($"Variance={dist2.Variance}");
            Console.WriteLine($"Skewness={dist2.Skewness}");
            Console.WriteLine($"Kurtosis={dist2.Kurtosis}");
            Console.WriteLine($"Entropy={dist2.Entropy}");
        }

        [TestMethod()]
        public void PDFTest() {
            Console.WriteLine(dist1);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist1.PDF(x)}");
            }
            Console.WriteLine(dist2);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist2.PDF(x)}");
            }
        }

        [TestMethod()]
        public void CDFTest() {
            //Console.WriteLine(dist1);
            //for (ddouble x = -4; x <= 4; x += 0.125) {
            //    Console.WriteLine($"cdf({x})={dist1.CDF(x)}");
            //}
            //Console.WriteLine(dist2);
            //for (ddouble x = -4; x <= 4; x += 0.125) {
            //    Console.WriteLine($"cdf({x})={dist2.CDF(x)}");
            //}
        }

        [TestMethod()]
        public void QuantileTest() {
            //Console.WriteLine(dist1);
            //for (int i = 0; i <= 10; i++) {
            //    ddouble p = (ddouble)i / 10;
            //    ddouble x = dist1.Quantile(p);

            //    Console.WriteLine($"quantile({p})={x}, cdf({x})={dist1.CDF(x)}");
            //}
            //Console.WriteLine(dist2);
            //for (int i = 0; i <= 10; i++) {
            //    ddouble p = (ddouble)i / 10;
            //    ddouble x = dist2.Quantile(p);

            //    Console.WriteLine($"quantile({p})={x}, cdf({x})={dist2.CDF(x)}");
            //}
        }
    }
}