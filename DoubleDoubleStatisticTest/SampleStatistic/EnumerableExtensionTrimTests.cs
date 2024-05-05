using DoubleDouble;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.SampleStatistic {
    [TestClass()]
    public class EnumerableExtensionTrimTests {
        readonly int[] testcase = [4, 8, 0, 7, 3, 10, 5, 9, 6, 1, 2];

        [TestMethod()]
        public void TrimTest() {
            CollectionAssert.AreEqual((new int[] { 4, 8, 7, 3, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).Trim(0.05, 0.95).ToArray());
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).Trim(0, 0.95).ToArray());
            CollectionAssert.AreEqual((new int[] { 4, 8, 7, 3, 10, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).Trim(0.05, 1).ToArray());
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 10, 5, 9, 6, 1, 2 }).Select(v => (double)v).ToArray(), testcase.Select(v => (double)v).Trim(0, 1).ToArray());

            CollectionAssert.AreEqual((new int[] { 4, 8, 7, 3, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).Trim(0.05, 0.95).ToArray());
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).Trim(0, 0.95).ToArray());
            CollectionAssert.AreEqual((new int[] { 4, 8, 7, 3, 10, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).Trim(0.05, 1).ToArray());
            CollectionAssert.AreEqual((new int[] { 4, 8, 0, 7, 3, 10, 5, 9, 6, 1, 2 }).Select(v => (ddouble)v).ToArray(), testcase.Select(v => (ddouble)v).Trim(0, 1).ToArray());
        }
    }
}