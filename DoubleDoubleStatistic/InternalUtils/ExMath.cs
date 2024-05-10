using DoubleDouble;
using System.Collections.ObjectModel;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.InternalUtils {
    internal static class ExMath {

        public static ddouble Pow3d2(ddouble x) {
            if (IsZero(x) || !IsFinite(x) || int.Abs(double.ILogB((double)x)) < 320) {
                return Sqrt(Cube(x));
            }
            else {
                return Cube(Sqrt(x));
            }
        }

        public static ddouble Pow2d3(ddouble x) {
            if (IsZero(x) || !IsFinite(x) || int.Abs(double.ILogB((double)x)) < 480) {
                return Cbrt(Square(x));
            }
            else {
                return Square(Cbrt(x));
            }
        }

        public static ReadOnlyCollection<ddouble> GenerateBesselIIntegerNuTable(ddouble x) {
            if (!IsFinite(x) || x > 256d) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            if (double.ILogB((double)x) < -240) {
                ddouble x2 = x * x;

                return new ReadOnlyCollection<ddouble>([
                    1d,
                    x * (8d + x2) / 16d,
                    x2 * (12d + x2) / 96d,
                    x * x2 / 48d
                ]);
            }

            const int n = 1024;

            ddouble[] vs = new ddouble[n + 2];
            vs[^1] = vs[^2] = double.ScaleB(1, -256);

            int truncate_n = n;

            for (int k = n - 1, m = 2 * n; k >= 0; k--, m -= 2) {
                vs[k] = m * vs[k + 1] / x + vs[k + 2];

                if (double.ILogB((double)vs[k]) > 256) {
                    for (int i = k; i < vs.Length; i++) {
                        vs[i] = Ldexp(vs[i], -256);

                        if (double.ILogB((double)vs[i]) <= -968) {
                            truncate_n = i;
                            break;
                        }
                    }
                }
            }

            vs = vs[..truncate_n];

            ddouble b0 = BesselI(0, x), r = b0 / vs[0];

            vs[0] = b0;

            for (int i = 1; i < vs.Length; i++) {
                vs[i] *= r;

                if (double.ILogB((double)vs[i]) <= -968) {
                    truncate_n = i;
                    break;
                }
            }

            vs = vs[..int.Min(n, truncate_n)];

            return new ReadOnlyCollection<ddouble>(vs);
        }
    }
}