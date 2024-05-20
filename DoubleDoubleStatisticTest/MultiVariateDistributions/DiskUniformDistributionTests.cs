using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.MultiVariateDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class DiskUniformDistributionTests {
        readonly DiskUniformDistribution dist = new();

        DiskUniformDistribution[] Dists => [
            dist,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (DiskUniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Covariance={dist.Covariance}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            Assert.IsTrue((new Vector(0, 0) - dist.Mean).Norm < 1e-5);
        }

        [TestMethod()]
        public void CovarianceTest() {
            Assert.IsTrue((new Matrix(new double[,]
                {{ 1,0}, {0,1}})
                - dist.Covariance).Norm < 1e-5);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.1447298858494002, (double)dist.Entropy, 1e-10);
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (DiskUniformDistribution dist in Dists) {
                double[][] samples = dist.Sample(random, count: 100000).ToArray();

                if (samples.Any(v => v[0] * v[0] + v[1] * v[1] - 1 > 1e-14)) {
                    Assert.Fail();
                }
            }
        }
    }
}
