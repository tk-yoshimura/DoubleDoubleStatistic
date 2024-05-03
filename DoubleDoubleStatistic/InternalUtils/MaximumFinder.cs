using DoubleDouble;

namespace DoubleDoubleStatistic.InternalUtils {
    internal static class MaximumFinder {
        public static ddouble BisectionFind(ddouble a, ddouble b, Func<ddouble, ddouble> f) {
            if (!(ddouble.IsFinite(a) && ddouble.IsFinite(b))) {
                return ddouble.NaN;
            }

            ddouble va = f(a), vb = f(b), h = b - a;

            for (int i = 0; i < 1024 && (ddouble.Abs(h) > ddouble.Abs(a) * 1e-30); i++) {
                ddouble xc = (a + b) * 0.5, vc = f(xc);

                if (va == vc && vb == vc) {
                    break;
                }
                if (ddouble.IsNaN(vc)) {
                    return ddouble.NaN;
                }

                ddouble xl = (a * 3d + b) * 0.25, vl = f(xl);
                ddouble xr = (a + b * 3d) * 0.25, vr = f(xr);

                int index = new (int index, ddouble v)[] { (0, va), (1, vl), (2, vc), (3, vr), (4, vb) }.MaxBy(item => item.v).index;

                if (index == 0) {
                    (b, vb) = (xl, vl);
                    h *= 0.25;
                }
                else if (index == 1) {
                    (b, vb) = (xc, vc);
                    h *= 0.5;
                }
                else if (index == 2) {
                    (a, va) = (xl, vl);
                    (b, vb) = (xr, vr);
                    h *= 0.5;
                }
                else if (index == 3) {
                    (a, va) = (xc, vc);
                    h *= 0.5;
                }
                else {
                    (a, va) = (xr, vr);
                    h *= 0.25;
                }
            }

            ddouble x = (a + b) * 0.5;

            return x;
        }
    }
}
