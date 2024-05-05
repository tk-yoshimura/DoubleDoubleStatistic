using DoubleDouble;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.SampleStatistic {
    [TestClass()]
    public class EnumerableExtensionMedianTests {
        readonly int[] testcase_1 = [4, 8, 0, 7, 3, 10, 5, 9, 6, 1, 2];
        readonly int[] testcase_2 = [4, 8, 0, 7, 3, 5, 9, 6, 1, 2];
        readonly int[] testcase_3 = [4];
        readonly int[] testcase_4 = [4, 8];

        [TestMethod()]
        public void MedianTest() {
            Assert.AreEqual(5, testcase_1.Select(v => (double)v).Median());
            Assert.AreEqual(4.5, testcase_2.Select(v => (double)v).Median());
            Assert.AreEqual(4, testcase_3.Select(v => (double)v).Median());
            Assert.AreEqual(6, testcase_4.Select(v => (double)v).Median());

            Assert.AreEqual(5, testcase_1.Select(v => (ddouble)v).Median());
            Assert.AreEqual(4.5, testcase_2.Select(v => (ddouble)v).Median());
            Assert.AreEqual(4, testcase_3.Select(v => (ddouble)v).Median());
            Assert.AreEqual(6, testcase_4.Select(v => (ddouble)v).Median());
        }
    }
}