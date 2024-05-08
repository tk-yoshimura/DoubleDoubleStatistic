namespace DoubleDoubleStatistic.RandomGeneration {
    public static partial class RandomExtension {
        private const UInt64 mask = 0x001F_FFFF_FFFF_FFFFuL;

        internal static double ToDouble(UInt64 n) {
            return double.ScaleB(n, -53);
        }

        internal static UInt64 NextBit53(this Random random) {
            UInt64 n = unchecked((UInt64)random.NextInt64()) & mask;

            return n;
        }

        public static double NextUniform(this Random random) {
            return double.ScaleB(unchecked((UInt64)random.NextInt64()), -63);
        }

        public static double NextUniformOpenInterval1(this Random random) {
            return ToDouble(NextBit53(random));
        }

        public static double NextUniformOpenInterval0(this Random random) {
            return ToDouble(NextBit53(random) + 1uL);
        }

        public static double NextUniformOpenInterval01(this Random random) {
            return ToDouble(NextBit53(random) | 1uL);
        }

        public static double NextGaussian(this Random random) {
            double x = random.NextUniformOpenInterval01(), y = random.NextUniformOpenInterval01();

            double z = double.Sqrt(-2d * double.Log(x)) * double.CosPi(2d * y);

            return z;
        }

        public static (double z1, double z2) NextGaussianX2(this Random random) {
            double x = random.NextUniformOpenInterval01(), y = random.NextUniformOpenInterval01();

            double lnx_sqrt = double.Sqrt(-2d * double.Log(x));
            (double s, double c) = double.SinCosPi(2d * y);

            double z1 = lnx_sqrt * s, z2 = lnx_sqrt * c;

            return (z1, z2);
        }
    }
}
