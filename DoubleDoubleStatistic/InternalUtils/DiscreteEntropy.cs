using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;

namespace DoubleDoubleStatistic.InternalUtils {
    internal static class DiscreteEntropy {
        public static ddouble Sum(DiscreteDistribution dist, int min, int max) {
            ddouble s = 0d, c = 1d;

            if (min <= max) {
                for (int k = min; k <= max; k++) {
                    ddouble p = dist.PMF(k);

                    if (p > 0d) {
                        ddouble ds = p * ddouble.Log(p);

                        s -= ds;

                        if (ddouble.Abs(ds) < ddouble.Abs(s) * 1e-30) {
                            break;
                        }
                    }

                    c -= p;

                    if (!(c >= 1e-29)) {
                        break;
                    }
                }
            }
            else {
                for (int k = min; k >= max; k--) {
                    ddouble p = dist.PMF(k);

                    if (p > 0d) {
                        ddouble ds = p * ddouble.Log(p);

                        s -= ds;

                        if (ddouble.Abs(ds) < ddouble.Abs(s) * 1e-30) {
                            break;
                        }
                    }

                    c -= p;

                    if (!(c >= 1e-29)) {
                        break;
                    }
                }
            }

            return s;
        }
    }
}
