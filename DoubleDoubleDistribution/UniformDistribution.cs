using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class UniformDistribution : Distribution {

        public ddouble Min { get; }
        public ddouble Max { get; }
        public ddouble Range { get; }

        private readonly ddouble pdf_norm;

        public UniformDistribution(ddouble min, ddouble max) {
            ValidateLocation(min);
            ValidateLocation(max);
            
            if (!(min < max)) {
                throw new ArgumentException($"Invalid location parameter. {min} < {max}");
            }

            this.Min = min;
            this.Max = max;
            this.Range = max - min;

            this.pdf_norm = 1d / this.Range;
        }

        public override ddouble PDF(ddouble x) {
            if (x < Min || x > Max) {
                return 0d;
            }

            ddouble pdf = this.pdf_norm;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < Min) {
                    return 0d;
                }

                if (x > Max) {
                    return 1d;
                }

                ddouble cdf = (x - Min) * this.pdf_norm;

                return cdf;
            }
            else {
                if (x < Min) {
                    return 1d;
                }

                if (x > Max) {
                    return 0d;
                }

                ddouble cdf = (Max - x) * this.pdf_norm;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble quantile = Min + p * this.Range;

                return quantile;
            }
            else {
                ddouble quantile = Max - p * this.Range;

                return quantile;
            }
        }

        public override (ddouble min, ddouble max) Support => (Min, Max);

        public override ddouble Mean => Ldexp(Min + Max, -1);
        public override ddouble Median => Ldexp(Min + Max, -1);
        
        public override ddouble Variance => Square(Range) / 12;
        public override ddouble Skewness => 0d;
        public override ddouble Kurtosis => -(ddouble)6 / 5;

        public override ddouble Entropy => Log(Range);

        public override string ToString() {
            return $"{typeof(UniformDistribution).Name}[min={Min},max={Max}]";
        }
    }
}
