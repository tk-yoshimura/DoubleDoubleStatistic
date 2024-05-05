using DoubleDouble;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.SampleStatistic {
    [TestClass()]
    public class EnumerableExtensionMomentTests {
        readonly int[] testcase = [22, 24, 26, 9, 33, 10, 23, 14, 19, 12, 15, 16];

        [TestMethod()]
        public void MeanTest() {
            Assert.AreEqual(18.583333333333332, testcase.Select(v => (double)v).Mean(), 1e-10);
            Assert.AreEqual(18.583333333333332, (double)testcase.Select(v => (ddouble)v).Mean(), 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(47.74305555555555, testcase.Select(v => (double)v).Variance(), 1e-10);
            Assert.AreEqual(47.74305555555555, (double)testcase.Select(v => (ddouble)v).Variance(), 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0.4375269495339061, testcase.Select(v => (double)v).Skewness(), 1e-10);
            Assert.AreEqual(0.4375269495339061, (double)testcase.Select(v => (ddouble)v).Skewness(), 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-0.6604750492561982, testcase.Select(v => (double)v).Kurtosis(), 1e-10);
            Assert.AreEqual(-0.6604750492561982, (double)testcase.Select(v => (ddouble)v).Kurtosis(), 1e-10);
        }

        [TestMethod()]
        public void StandardDeviationTest() {
            Assert.AreEqual(6.909634979907082, testcase.Select(v => (double)v).StandardDeviation(), 1e-10);
            Assert.AreEqual(6.909634979907082, (double)testcase.Select(v => (ddouble)v).StandardDeviation(), 1e-10);
        }
    }
}