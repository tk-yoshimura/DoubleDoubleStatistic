using DoubleDouble;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.SampleStatistic {
    [TestClass()]
    public class EnumerableExtensionIndexOfTests {
        readonly int[] testcase = [4, 8, 0, 7, 3, 10, 5, 9, 6, 1, 2];

        [TestMethod()]
        public void MedianTest() {
            Assert.AreEqual(5, testcase.Select(v => (double)v).MaxIndexOf());
            Assert.AreEqual(2, testcase.Select(v => (double)v).MinIndexOf());

            Assert.AreEqual(5, testcase.Select(v => (ddouble)v).MaxIndexOf());
            Assert.AreEqual(2, testcase.Select(v => (ddouble)v).MinIndexOf());
        }
    }
}