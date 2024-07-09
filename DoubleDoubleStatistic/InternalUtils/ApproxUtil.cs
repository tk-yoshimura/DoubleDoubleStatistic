using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic.InternalUtils {
    internal static class ApproxUtil {

        public static ddouble Pade(ddouble x, ReadOnlyCollection<(ddouble c, ddouble d)> table) {
            (ddouble sc, ddouble sd) = table[^1];

            for (int i = table.Count - 2; i >= 0; i--) {
                sc = sc * x + table[i].c;
                sd = sd * x + table[i].d;
            }

            Debug.Assert(sd >= 0.5, $"pade denom digits loss! {x}");

            return sc / sd;
        }

        public static ddouble Poly(ddouble x, ReadOnlyCollection<ddouble> table) {
            ddouble s = table[^1];

            for (int i = table.Count - 2; i >= 0; i--) {
                s = s * x + table[i];
            }

            return s;
        }

        public static double Pade(double x, (ReadOnlyCollection<double> c, ReadOnlyCollection<double> d) table) {
            (double sc, double sd) = (Poly(x, table.c), Poly(x, table.d));

            Debug.Assert(sd >= 0.5, $"pade denom digits loss! {x}");

            return sc / sd;
        }

        public static double Poly(double x, ReadOnlyCollection<double> table) {
            double s = table[^1];

            for (int i = table.Count - 2; i >= 0; i--) {
                s = s * x + table[i];
            }

            return s;
        }
    }
}