using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.Utils {
    [TestClass()]
    public class QuantileLinearFitterTests {
        [TestMethod()]
        public void FitTest1() {
            NormalDistribution dist = new(2, 7);

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            NormalDistribution dist_fit = QuantileLinearFitter<NormalDistribution>.Fit(new NormalDistribution(), samples, (0.125, 0.875));

            Console.WriteLine(dist_fit);

            Assert.AreEqual((double)dist.Mu, (double)dist_fit.Mu, 0.01);
            Assert.AreEqual((double)dist.Sigma, (double)dist_fit.Sigma, 0.01);
        }

        [TestMethod()]
        public void FitTest2() {
            LandauDistribution dist = new(2, 7);

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            LandauDistribution dist_fit = QuantileLinearFitter<LandauDistribution>.Fit(new LandauDistribution(), samples, (0.125, 0.875));

            Console.WriteLine(dist_fit);

            Assert.AreEqual((double)dist.Mu, (double)dist_fit.Mu, 0.01);
            Assert.AreEqual((double)dist.C, (double)dist_fit.C, 0.01);
        }
    }
}