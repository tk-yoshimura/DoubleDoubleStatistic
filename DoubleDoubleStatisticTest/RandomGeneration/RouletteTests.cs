using DoubleDoubleStatistic.RandomGeneration;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.RandomGeneration {
    [TestClass()]
    public class RouletteTests {
        [TestMethod()]
        public void RouletteTest() {
            Random random = new(1234);

            double[] probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

            for (int n = 1; n < probs.Length; n++) {
                Roulette roulette = new(probs[..n]);

                int[] samples = roulette.NextIndex(random, count: 100000).ToArray();

                for (int i = 0; i < n; i++) {
                    Assert.AreEqual(probs[i] / probs[..n].Sum(), samples.Count(c => c == i) / (double)samples.Length, 0.01, $"{n},{i}");
                }
            }
        }
    }
}