using DoubleDouble;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.SampleStatistic {
    [TestClass()]
    public class EnumerableExtensionQuantileTests {
        readonly int[] testcase_1 = [4, 8, 0, 7, 3, 10, 5, 9, 6, 1, 2];
        readonly int[] testcase_2 = [4];

        [TestMethod()]
        public void QuantileTest() {
            for (int i = 0; i <= 40; i++) {
                double p = i / 40d;

                Assert.AreEqual(p * 10, testcase_1.Select(v => (double)v).Quantile(p), 1e-10, $"{p}");
                Assert.AreEqual(p * 10, (double)testcase_1.Select(v => (ddouble)v).Quantile(p), 1e-10, $"{p}");
            }

            for (int i = 0; i <= 40; i++) {
                double p = i / 40d;

                Assert.AreEqual(4, testcase_2.Select(v => (double)v).Quantile(p), 1e-10, $"{p}");
                Assert.AreEqual(4, (double)testcase_2.Select(v => (ddouble)v).Quantile(p), 1e-10, $"{p}");
            }
        }

        [TestMethod()]
        public void QuantileTupleTest() {
            Assert.AreEqual((2.5, 7.5), testcase_1.Select(v => (double)v).Quantile(0.25, 0.75));
            Assert.AreEqual((2.5, 7.5), testcase_1.Select(v => (ddouble)v).Quantile(0.25, 0.75));

            Assert.AreEqual((2.5, 5.0, 7.5), testcase_1.Select(v => (double)v).Quantile(0.25, 0.50, 0.75));
            Assert.AreEqual((2.5, 5.0, 7.5), testcase_1.Select(v => (ddouble)v).Quantile(0.25, 0.50, 0.75));

            Assert.AreEqual((1.25, 2.5, 5.0, 7.5), testcase_1.Select(v => (double)v).Quantile(0.125, 0.25, 0.50, 0.75));
            Assert.AreEqual((1.25, 2.5, 5.0, 7.5), testcase_1.Select(v => (ddouble)v).Quantile(0.125, 0.25, 0.50, 0.75));
        }
    }
}