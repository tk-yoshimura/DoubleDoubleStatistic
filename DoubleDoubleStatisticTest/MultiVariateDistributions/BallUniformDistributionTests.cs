using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.MultiVariateDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class BallUniformDistributionTests {
        readonly BallUniformDistribution dist = new();

        BallUniformDistribution[] Dists => [
            dist,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BallUniformDistribution dist in Dists) {
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
            Assert.IsTrue((new Vector(0, 0, 0) - dist.Mean).Norm < 1e-5);
        }

        [TestMethod()]
        public void CovarianceTest() {
            Assert.IsTrue((new Matrix(new double[,]
                {{ 1,0,0}, {0,1,0}, {0, 0, 1}}) / 5
                - dist.Covariance).Norm < 1e-5);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.4324119583011812, (double)dist.Entropy, 1e-10);
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (BallUniformDistribution dist in Dists) {
                double[][] samples = dist.Sample(random, count: 100000).ToArray();

                (Vector mu, Matrix matrix) = samples.Covariance();

                if (samples.Any(v => v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - 1 > 1e-14)) {
                    Assert.Fail();
                }
            }
        }
    }
}
