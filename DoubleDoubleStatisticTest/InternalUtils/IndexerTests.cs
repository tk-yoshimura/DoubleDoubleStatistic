using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.ObjectModel;

namespace DoubleDoubleStatisticTest.InternalUtils {
    [TestClass()]
    public class IndexerTests {
        [TestMethod()]
        public void BisectionSearchTest1() {
            ReadOnlyCollection<ddouble> thresholds = new([0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9]);

            Assert.AreEqual(0, Indexer.BisectionSearch(0.15, thresholds));
            Assert.AreEqual(5, Indexer.BisectionSearch(0.85, thresholds));
            Assert.AreEqual(1, Indexer.BisectionSearch(0.2, thresholds));
            Assert.AreEqual(1, Indexer.BisectionSearch(0.3, thresholds));
            Assert.AreEqual(2, Indexer.BisectionSearch(0.4, thresholds));
            Assert.AreEqual(3, Indexer.BisectionSearch(0.6, thresholds));
            Assert.AreEqual(0, Indexer.BisectionSearch(0.1, thresholds));
            Assert.AreEqual(5, Indexer.BisectionSearch(0.9, thresholds));

            Assert.AreEqual(0, Indexer.BisectionSearch(0.0, thresholds));
            Assert.AreEqual(5, Indexer.BisectionSearch(1.0, thresholds));
        }

        [TestMethod()]
        public void BisectionSearchTest2() {
            ReadOnlyCollection<double> thresholds = new([0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9]);

            Assert.AreEqual(0, Indexer.BisectionSearch(0.15, thresholds));
            Assert.AreEqual(5, Indexer.BisectionSearch(0.85, thresholds));
            Assert.AreEqual(1, Indexer.BisectionSearch(0.2, thresholds));
            Assert.AreEqual(1, Indexer.BisectionSearch(0.3, thresholds));
            Assert.AreEqual(2, Indexer.BisectionSearch(0.4, thresholds));
            Assert.AreEqual(3, Indexer.BisectionSearch(0.6, thresholds));
            Assert.AreEqual(0, Indexer.BisectionSearch(0.1, thresholds));
            Assert.AreEqual(5, Indexer.BisectionSearch(0.9, thresholds));

            Assert.AreEqual(0, Indexer.BisectionSearch(0.0, thresholds));
            Assert.AreEqual(5, Indexer.BisectionSearch(1.0, thresholds));
        }
    }
}
