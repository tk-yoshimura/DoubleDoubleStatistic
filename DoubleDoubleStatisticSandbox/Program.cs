using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;

namespace DoubleDoubleStatisticSandbox {
    internal class Program {
        static void Main() {
            QGaussianDistribution dist = new(2.5, 2);

            for (ddouble x = -16; x <= 16; x += 0.0625) {
                Console.WriteLine($"{x},{dist.PDF(x)}");
            }

            dist.Quantile("0.4");

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
