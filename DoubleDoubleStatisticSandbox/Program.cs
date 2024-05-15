using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;

namespace DoubleDoubleStatisticSandbox {
    internal class Program {
        static void Main() {
            Random random = new(1234);
            NakagamiDistribution dist = new(3, 4);

            double[] samples = dist.Sample(random, 10000).ToArray();

            NakagamiDistribution? dist_fit = NakagamiDistribution.Fit(samples, (0.01, 0.99)).dist;

            for (ddouble x = -16; x <= 16; x += 0.0625) {
                Console.WriteLine($"{x},{dist.PDF(x)},{dist_fit?.PDF(x)}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
