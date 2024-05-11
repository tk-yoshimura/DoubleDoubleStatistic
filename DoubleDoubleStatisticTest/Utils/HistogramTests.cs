using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.Utils {
    [TestClass()]
    public class HistogramTests {
        [TestMethod()]
        public void HistogramTest1() {
            NormalDistribution dist = new();

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            Histogram hist = new(samples, 32, (-4, 4));

            for (int i = 0; i <= hist.Bins; i++) {
                Assert.AreEqual(-4 + i * 0.25, hist.BinsEdge[i]);
            }

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual(-3.875 + i * 0.25, hist.BinsCentroid[i]);
            }

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentroid[i]), (double)hist.Density[i], 0.01);
            }
        }

        [TestMethod()]
        public void HistogramTest2() {
            ChiSquareDistribution dist = new(3);

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            Histogram hist = new(samples, 64, (1, 9));

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentroid[i]), (double)hist.Density[i], 0.01);
            }
        }
    }
}