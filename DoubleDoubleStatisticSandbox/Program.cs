using DoubleDoubleStatistic.ContinuousDistributions;

namespace DoubleDoubleStatisticSandbox {
    internal class Program {
        static void Main() {
            Random random = new();
            VoigtDistribution dist = new(1, 1);

            for (int i = 0; i < 1024; i++) {
                double y = dist.Sample(random);

                Console.WriteLine($"{y}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
