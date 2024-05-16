using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.MultiVariateDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class MultiVariateNormalDistribution1DTests {
        readonly MultiVariateNormalDistribution dist = new(
            mu: new Vector(0.5),
            cov: new Matrix(new double[,] { { 4.0 } })
        );

        readonly NormalDistribution dist_expected = new NormalDistribution(0.5, 2);

        [TestMethod()]
        public void InfoTest() {
            Console.WriteLine(dist);
            Console.WriteLine($"Mu={dist.Mu}");
            Console.WriteLine($"Mean={dist.Mean}");
            Console.WriteLine($"Mode={dist.Mode}");
            Console.WriteLine($"Covariance={dist.Covariance}");
            Console.WriteLine($"Entropy={dist.Entropy}");
            Console.WriteLine(dist.Formula);
        }

        [TestMethod()]
        public void MeanTest() {
            Assert.AreEqual(dist_expected.Mean, dist.Mean[0]);
        }

        [TestMethod()]
        public void ModeTest() {
            Assert.AreEqual(dist_expected.Mode, dist.Mode[0]);
        }

        [TestMethod()]
        public void CovarianceTest() {
            Assert.AreEqual(dist_expected.Variance, dist.Covariance[0, 0]);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual((double)dist_expected.Entropy, (double)dist.Entropy);
        }

        [TestMethod()]
        public void PDFTest() {
            Console.WriteLine(dist);
            for (int x = -5; x <= 5; x++) {
                ddouble pdf = dist.PDF(new Vector(x));

                Console.WriteLine($"pdf({x})={pdf}");
            }
        }

        [TestMethod()]
        public void RandomGenerateAndFitTest() {
            Random random = new(1234);

            double[][] samples = dist.Sample(random, count: 100000).ToArray();

            MultiVariateNormalDistribution? dist_fit = MultiVariateNormalDistribution.Fit(samples);

            Assert.IsNotNull(dist_fit);

            Console.WriteLine(dist_fit);

            Assert.IsTrue((dist.Mu - dist_fit.Mu).Norm < 1e-2, $"{dist},mu");
            Assert.IsTrue((dist.Covariance - dist_fit.Covariance).Norm < 1e-1, $"{dist},cov");
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            for (int x = -5; x <= 5; x++) {
                ddouble expected = dist_expected.PDF(x);
                ddouble actual = dist.PDF(new Vector(x));

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
}
