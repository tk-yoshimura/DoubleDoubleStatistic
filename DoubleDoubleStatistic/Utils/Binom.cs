using System.Numerics;

namespace DoubleDoubleStatistic {
    internal static class Binom {
        static readonly Dictionary<(int n, int k), BigInteger> table = [];

        public static BigInteger Value(int n, int k) {
            if (n < 0 || k > n) {
                return BigInteger.Zero;
            }

            if (k == 0 || k == n) {
                return BigInteger.One;
            }

            if (table.TryGetValue((n, k), out BigInteger value)) {
                return value;
            }

            BigInteger new_value = Value(n - 1, k - 1) + Value(n - 1, k);

            table.Add((n, k), new_value);

            return new_value;
        }
    }
}
