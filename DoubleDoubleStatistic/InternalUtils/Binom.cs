using System.Collections.Concurrent;
using System.Numerics;

namespace DoubleDoubleStatistic.InternalUtils {
    internal static class Binom {
        static readonly ConcurrentDictionary<(int n, int k), BigInteger> table = [];

        public static BigInteger Value(int n, int k) {
            if (n < 0 || k > n) {
                return BigInteger.Zero;
            }

            if (k == 0 || k == n) {
                return BigInteger.One;
            }

            if (!table.TryGetValue((n, k), out BigInteger value)) {
                value = Value(n - 1, k - 1) + Value(n - 1, k);
                table[(n, k)] = value;
            }

            return value;
        }
    }
}
