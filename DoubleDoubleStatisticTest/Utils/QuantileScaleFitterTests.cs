using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.Utils {
    [TestClass()]
    public class QuantileScaleFitterTests {
        [TestMethod()]
        public void FitTest1() {
            NormalDistribution dist = new(sigma: 7);

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            NormalDistribution dist_fit = QuantileScaleFitter<NormalDistribution>.Fit(new NormalDistribution(), samples, (0.125, 0.875));

            Console.WriteLine(dist_fit);

            Assert.AreEqual(0d, dist_fit.Mu);
            Assert.AreEqual((double)dist.Sigma, (double)dist_fit.Sigma, 0.01);
        }

        [TestMethod()]
        public void FitTest2() {
            ExponentialDistribution dist = new(7);

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            ExponentialDistribution dist_fit = QuantileScaleFitter<ExponentialDistribution>.Fit(new ExponentialDistribution(), samples, (0.125, 0.875));

            Console.WriteLine(dist_fit);

            Assert.AreEqual((double)dist.Theta, (double)dist_fit.Theta, 0.01);
        }
    }
}