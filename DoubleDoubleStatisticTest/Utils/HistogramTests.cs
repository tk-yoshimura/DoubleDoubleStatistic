using DoubleDouble;
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
                Assert.AreEqual(-3.875 + i * 0.25, hist.BinsCentor[i]);
            }

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentor[i]), (double)hist.Density[i], 0.01);
            }
        }

        [TestMethod()]
        public void HistogramTest2() {
            ChiSquareDistribution dist = new(3);

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            Histogram hist = new(samples, 64, (1, 9));

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentor[i]), (double)hist.Density[i], 0.01);
            }
        }

        [TestMethod()]
        public void HistogramTest3() {
            NormalDistribution dist = new();

            double[] samples = dist.Sample(new Random(1234), 1000000).ToArray();

            Histogram hist = new(samples, new double[] { -1, 0, 2, 4, 8 });

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentor[i]), (double)hist.Density[i], 0.02);
            }
        }

        [TestMethod()]
        public void HistogramTest4() {
            NormalDistribution dist = new();

            ddouble[] samples = dist.Sample(new Random(1234), 1000000).Select(v => (ddouble)v).ToArray();

            Histogram hist = new(samples, 32, (-4, 4));

            for (int i = 0; i <= hist.Bins; i++) {
                Assert.AreEqual(-4 + i * 0.25, hist.BinsEdge[i]);
            }

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual(-3.875 + i * 0.25, hist.BinsCentor[i]);
            }

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentor[i]), (double)hist.Density[i], 0.01);
            }
        }

        [TestMethod()]
        public void HistogramTest5() {
            ChiSquareDistribution dist = new(3);

            ddouble[] samples = dist.Sample(new Random(1234), 1000000).Select(v => (ddouble)v).ToArray();

            Histogram hist = new(samples, 64, (1, 9));

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentor[i]), (double)hist.Density[i], 0.01);
            }
        }

        [TestMethod()]
        public void HistogramTest6() {
            NormalDistribution dist = new();

            ddouble[] samples = dist.Sample(new Random(1234), 1000000).Select(v => (ddouble)v).ToArray();

            Histogram hist = new(samples, new ddouble[] { -1, 0, 2, 4, 8 });

            for (int i = 0; i < hist.Bins; i++) {
                Assert.AreEqual((double)dist.PDF(hist.BinsCentor[i]), (double)hist.Density[i], 0.02);
            }
        }
    }
}