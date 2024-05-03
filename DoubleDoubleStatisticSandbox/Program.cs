using DoubleDouble;
using DoubleDoubleStatistic;

namespace DoubleDoubleStatisticSandbox {
    internal class Program {
        static void Main() {
            LandauDistribution dist = new(1);

            for (ddouble x = -8; x <= 24; x += 1 / 32d) {
                ddouble y = dist.PDF(x);

                Console.WriteLine($"{x},{y}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
