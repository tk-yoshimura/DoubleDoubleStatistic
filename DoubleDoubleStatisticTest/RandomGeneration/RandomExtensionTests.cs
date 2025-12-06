using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.RandomGeneration {
    [TestClass()]
    public class RandomExtensionTests {
        [TestMethod()]
        public void ToDoubleTest() {
            Assert.AreEqual(0d, RandomExtension.ToDouble(0));
            Assert.IsGreaterThan(0, RandomExtension.ToDouble(1));
            Assert.IsLessThan(1, RandomExtension.ToDouble(0x001F_FFFF_FFFF_FFFFuL));
            Assert.AreEqual(1d, RandomExtension.ToDouble(0x0020_0000_0000_0000uL));
            Assert.AreEqual(1d, double.BitIncrement(RandomExtension.ToDouble(0x001F_FFFF_FFFF_FFFFuL)));
            Assert.AreEqual(1d, RandomExtension.ToDouble(1) + RandomExtension.ToDouble(0x001F_FFFF_FFFF_FFFFuL));
        }

        [TestMethod()]
        public void NextBit53Test() {
            Random random = new(1234);

            for (int i = 0; i < 1024; i++) {
                UInt64 n = RandomExtension.NextBit53(random);

                Assert.IsGreaterThanOrEqualTo(0uL, n);
                Assert.IsLessThanOrEqualTo(0x001F_FFFF_FFFF_FFFFuL, n);
            }
        }

        [TestMethod()]
        public void NextUniformTest() {
            Random random = new(1234);

            double sv = 0d;

            for (int i = 0; i < 1024; i++) {
                double v = random.NextUniform();

                sv += v;

                Assert.IsGreaterThanOrEqualTo(0d, v);
                Assert.IsLessThanOrEqualTo(1d, v);
            }

            Assert.AreEqual(512, sv, 16);
        }

        [TestMethod()]
        public void NextUniformOpenInterval0Test() {
            Random random = new(1234);

            double sv = 0d;

            for (int i = 0; i < 1024; i++) {
                double v = random.NextUniformOpenInterval0();

                sv += v;

                Assert.IsGreaterThan(0d, v);
                Assert.IsLessThanOrEqualTo(1d, v);
            }

            Assert.AreEqual(512, sv, 16);
        }

        [TestMethod()]
        public void NextUniformOpenInterval1Test() {
            Random random = new(1234);

            double sv = 0d;

            for (int i = 0; i < 1024; i++) {
                double v = random.NextUniformOpenInterval1();

                sv += v;

                Assert.IsGreaterThanOrEqualTo(0d, v);
                Assert.IsLessThan(1d, v);
            }

            Assert.AreEqual(512, sv, 16);
        }

        [TestMethod()]
        public void NextUniformOpenInterval01Test() {
            Random random = new(1234);

            double sv = 0d;

            for (int i = 0; i < 1024; i++) {
                double v = random.NextUniformOpenInterval01();

                sv += v;

                Assert.IsGreaterThan(0d, v);
                Assert.IsLessThan(1d, v);
            }

            Assert.AreEqual(512, sv, 16);
        }

        [TestMethod()]
        public void NextGaussianTest() {
            Random random = new(1234);

            List<double> vs = [];

            for (int i = 0; i < 1024; i++) {
                double v = random.NextGaussian();

                vs.Add(v);
            }

            Assert.AreEqual(0d, vs.Mean(), 4);
            Assert.AreEqual(1d, vs.Variance(), 0.25);
            Assert.AreEqual(0d, vs.Skewness(), 0.25);
            Assert.AreEqual(0d, vs.Kurtosis(), 0.25);
        }

        [TestMethod()]
        public void NextGaussianX2Test() {
            Random random = new(1234);

            List<double> vs = [];

            for (int i = 0; i < 1024; i++) {
                (double v1, double v2) = random.NextGaussianX2();

                vs.Add(v1);
                vs.Add(v2);
            }

            Assert.AreEqual(0d, vs.Mean(), 4);
            Assert.AreEqual(1d, vs.Variance(), 0.25);
            Assert.AreEqual(0d, vs.Skewness(), 0.25);
            Assert.AreEqual(0d, vs.Kurtosis(), 0.25);
        }
    }
}