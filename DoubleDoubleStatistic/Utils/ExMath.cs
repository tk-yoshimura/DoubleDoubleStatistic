using DoubleDouble;

namespace DoubleDoubleStatistic {
    internal static class ExMath {

        public static ddouble Pow3d2(ddouble x) {
            if (ddouble.IsZero(x) || !ddouble.IsFinite(x) || int.Abs(double.ILogB((double)x)) < 320) {
                return ddouble.Sqrt(ddouble.Cube(x));
            }
            else {
                return ddouble.Cube(ddouble.Sqrt(x));
            }
        }

        public static ddouble Pow2d3(ddouble x) {
            if (ddouble.IsZero(x) || !ddouble.IsFinite(x) || int.Abs(double.ILogB((double)x)) < 480) {
                return ddouble.Cbrt(ddouble.Square(x));
            }
            else {
                return ddouble.Square(ddouble.Cbrt(x));
            }
        }
    }
}