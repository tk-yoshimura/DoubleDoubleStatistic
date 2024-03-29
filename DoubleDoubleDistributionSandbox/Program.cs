using DoubleDouble;
using DoubleDoubleDistribution;

namespace DoubleDoubleDistributionSandbox {
    internal class Program {
        static void Main() {
            LandauDistribution dist = new(mu: 0, c: ddouble.PI / 2);

            for (ddouble x = -4; x <= 16; x += 1 / 256d) {
                ddouble y = dist.PDF(x);

                Console.WriteLine($"{x},{y}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
