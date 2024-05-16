using Algebra;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.MultiVariateDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class MultiVariateNormalDistribution4DTests {
        [TestMethod()]
        public void RandomGenerateAndFitTest() {
            Random random = new(1234);

            UniformDistribution uniform = new(-1, 1);

            for (int i = 0; i < 16; i++) {
                Vector mu = new Vector(uniform.Sample(random, 4).ToArray());
                Matrix m = Matrix.HConcat(
                    uniform.Sample(random, 4).ToArray(),
                    uniform.Sample(random, 4).ToArray(),
                    uniform.Sample(random, 4).ToArray(),
                    uniform.Sample(random, 4).ToArray()
                );

                Matrix cov = m * m.T;

                MultiVariateNormalDistribution dist = new(mu, cov);

                double[][] samples = dist.Sample(random, count: 100000).ToArray();

                MultiVariateNormalDistribution? dist_fit = MultiVariateNormalDistribution.Fit(samples);

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist);
                Console.WriteLine(dist_fit);

                Assert.IsTrue((dist.Mu - dist_fit.Mu).Norm < 1e-1, $"{dist},mu");
                Assert.IsTrue((dist.Covariance - dist_fit.Covariance).Norm < 1e-1, $"{dist},cov");
            }
        }
    }
}
