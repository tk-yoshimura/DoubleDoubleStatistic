using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class WignerSemicircleDistribution : Distribution {

        public ddouble Radius { get; }

        private readonly ddouble pdf_norm, radius_sq;

        public WignerSemicircleDistribution(ddouble radius) {
            ValidateScale(radius);

            this.Radius = radius;

            this.radius_sq = radius * radius;
            this.pdf_norm = 1d / (PI * radius_sq);
        }

        public override ddouble PDF(ddouble x) {
            if (x < -Radius || x > Radius) {
                return 0d;
            }

            ddouble pdf = 2 * this.pdf_norm * Sqrt(radius_sq - x * x);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < -Radius) {
                    return 0d;
                }

                if (x > Radius) {
                    return 1d;
                }

                ddouble cdf = 0.5d + x * this.pdf_norm * Sqrt(radius_sq - x * x) + Asin(x / Radius) * RcpPI;

                return cdf;
            }
            else {
                if (x < -Radius) {
                    return 1d;
                }

                if (x > Radius) {
                    return 0d;
                }

                ddouble cdf = 0.5d - x * this.pdf_norm * Sqrt(radius_sq - x * x) - Asin(x / Radius) * RcpPI;

                return cdf;
            }
        }

        public override (ddouble min, ddouble max) Support => (-Radius, Radius);

        public override ddouble Mean => 0d;
        public override ddouble Median => 0d;

        public override ddouble Variance => radius_sq / 4;
        public override ddouble Skewness => 0d;
        public override ddouble Kurtosis => -1;

        public override ddouble Entropy => Log(PI * Radius) - 0.5d;

        public override string ToString() {
            return $"{typeof(WignerSemicircleDistribution).Name}[radius={Radius}]";
        }
    }
}
