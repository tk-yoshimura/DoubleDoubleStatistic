using DoubleDouble;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.SampleStatistic {
    [TestClass()]
    public class EnumerableExtensionNaNTests {
        readonly double[] testcase = [4, 8, 0, 7, 3, double.NaN, 5, 9, 6, 1, 2];

        [TestMethod()]
        public void NaNAsTest() {
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, -1, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).NaNAs(-1).ToArray());

            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, -1, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).NaNAs(-1).ToArray());
        }

        [TestMethod()]
        public void NaNAsZeroTest() {
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 0, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).NaNAsZero().ToArray());

            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 0, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).NaNAsZero().ToArray());
        }

        [TestMethod()]
        public void TrimNaNTest() {
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).TrimNaN().ToArray());

            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).TrimNaN().ToArray());
        }
    }
}