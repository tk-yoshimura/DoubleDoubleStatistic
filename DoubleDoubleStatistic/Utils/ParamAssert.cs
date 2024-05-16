using DoubleDouble;

namespace DoubleDoubleStatistic.Utils {
    public static class ParamAssert {

        public static bool IsFinitePositive(ddouble x) => x > 0d && ddouble.IsFinite(x);

        public static void ValidateSupport(string param_name, bool condition) {
            if (!condition) {
                throw new ArgumentOutOfRangeException(param_name, "Invalid support parameter.");
            }
        }

        public static void ValidateScale(string param_name, bool condition) {
            if (!condition) {
                throw new ArgumentOutOfRangeException(param_name, "Invalid scale parameter.");
            }
        }

        public static void ValidateLocation(string param_name, bool condition) {
            if (!condition) {
                throw new ArgumentOutOfRangeException(param_name, "Invalid location parameter.");
            }
        }

        public static void ValidateShape(string param_name, bool condition) {
            if (!condition) {
                throw new ArgumentOutOfRangeException(param_name, "Invalid shape parameter.");
            }
        }

        public static void ValidateNonCentricity(string param_name, bool condition) {
            if (!condition) {
                throw new ArgumentOutOfRangeException(param_name, "Invalid non-centricity parameter.");
            }
        }
    }
}
