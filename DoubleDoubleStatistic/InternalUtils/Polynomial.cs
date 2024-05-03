using System.Collections.ObjectModel;
using System.Numerics;

namespace DoubleDoubleStatistic.InternalUtils {
    internal class Polynomial {
        public ReadOnlyCollection<BigInteger> Coef { get; }

        public int Degree => Coef.Count - 1;

        public Polynomial(BigInteger[] coef) {
            Coef = new ReadOnlyCollection<BigInteger>(coef);
        }

        public static Polynomial Zero => new([]);

        public BigInteger this[int index] => index < Coef.Count ? Coef[index] : BigInteger.Zero;

        public static Polynomial XMul(Polynomial p, int degree) {
            ArgumentOutOfRangeException.ThrowIfNegative(degree, nameof(degree));

            BigInteger[] coef = Enumerable.Repeat(BigInteger.Zero, degree).Concat(p.Coef).ToArray();

            return new(coef);
        }

        public static Polynomial XShift(Polynomial p, BigInteger c) {
            Polynomial q = Zero;
            Polynomial x = new([c, 1]), v = new([1]);

            for (int i = 0; i < p.Coef.Count; i++) {
                q += v * p.Coef[i];
                v *= x;
            }

            return q;
        }

        public static Polynomial operator *(BigInteger n, Polynomial p) {
            return new(p.Coef.Select(c => n * c).ToArray());
        }

        public static Polynomial operator *(Polynomial p, BigInteger n) => n * p;

        public static Polynomial operator *(Polynomial p1, Polynomial p2) {
            Polynomial p = Zero;

            for (int i = 0; i < p2.Coef.Count; i++) {
                p += XMul(p1 * p2.Coef[i], i);
            }

            return p;
        }

        public static Polynomial operator +(Polynomial p1, Polynomial p2) {
            BigInteger[] coef = new BigInteger[Math.Max(p1.Coef.Count, p2.Coef.Count)];

            for (int i = 0; i < coef.Length; i++) {
                if (i < p1.Coef.Count && i < p2.Coef.Count) {
                    coef[i] = p1.Coef[i] + p2.Coef[i];
                }
                else if (i < p1.Coef.Count) {
                    coef[i] = p1.Coef[i];
                }
                else if (i < p2.Coef.Count) {
                    coef[i] = p2.Coef[i];
                }
            }

            return new(coef);
        }

        public static Polynomial operator -(Polynomial p1, Polynomial p2) => p1 + (-p2);

        public static Polynomial operator +(Polynomial p) => p;

        public static Polynomial operator -(Polynomial p) {
            return new(p.Coef.Select(BigInteger.Negate).ToArray());
        }

        public override string ToString() {
            if (Degree < 0) {
                return "0";
            }

            List<string> strs = ["+"];

            if (Coef[0] != 0) {
                strs.Add($"{Coef[0]}");
            }

            if (Degree >= 1) {
                if (Coef[1] != 0) {
                    strs.Add($"{Coef[1]}x");
                }

                for (int i = 2; i < Coef.Count; i++) {
                    if (Coef[i] != 0) {
                        strs.Add($"{Coef[i]}x^{i}");
                    }
                }
            }

            string str = string.Join("+", strs).Replace("+1x", "+x").Replace("-1x", "-x").Replace("+-", "-").Replace("+-", "-").Replace("++", "");

            if (string.IsNullOrEmpty(str) || str == "+") {
                return "0";
            }

            return str;
        }
    }
}
