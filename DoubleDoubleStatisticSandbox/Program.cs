using DoubleDouble;
using DoubleDoubleStatistic;

namespace DoubleDoubleStatisticSandbox {
    internal class Program {
        static void Main() {
            BradfordDistribution dist = new(1);

            for (ddouble x = 0; x <= 1; x += 1 / 256d) {
                ddouble y = dist.PDF(x);

                Console.WriteLine($"{x},{y}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
