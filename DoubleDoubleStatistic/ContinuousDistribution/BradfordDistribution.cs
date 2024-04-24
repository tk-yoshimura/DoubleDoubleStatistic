using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class BradfordDistribution : ContinuousDistribution {

        public ddouble C { get; }

        private readonly ddouble pdf_norm;

        public BradfordDistribution(ddouble c) {
            C = c;
            pdf_norm = c / Log1p(c);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x) || x > 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm / (1d + C * x);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }

                if (x > 1d) {
                    return 1d;
                }

                ddouble cdf = Log1p(C * x) * pdf_norm / C;

                return cdf;
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                if (x > 1d) {
                    return 0d;
                }

                ddouble cdf = 1d - Log1p(C * x) * pdf_norm / C;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return Quantile(1d - p);
            }

            ddouble x = Expm1(p * C / pdf_norm) / C;

            x = Clamp(x, 0d, 1d);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (0d, 1d);

        public override ddouble Mean {
            get {
                ddouble k = Log1p(C);
                return (C - k) / (C * k);
            }
        }

        public override ddouble Median => Expm1(C / (pdf_norm * 2)) / C;

        public override ddouble Mode => 0d;

        public override ddouble Variance {
            get {
                ddouble k = Log1p(C);
                return ((C + 2d) * k - 2d * C) / (2d * C * k * k);
            }
        }

        public override ddouble Skewness {
            get {
                ddouble k = Log1p(C);
                return Sqrt2 * (2d * k * k * (3d + C * (3d + C)) - 9d * C * k * (C + 2d) + 12d * C * C) /
                    (Sqrt(C * (2 * k + C * (k - 2d))) * (6 * k + 3d * C * (k - 2d)));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble k = Log1p(C);
                return (12d * k * k * k
                    + 6d * C * k * k * (3d * k - 14d)
                    + 12d * C * C * k * (k - 4d) * (k - 3d)
                    + C * C * C * (k - 3d) * (3d * k * k - 16d * k + 24d))
                    / (3 * C * Square(2d * k + C * (k - 2d)));
            }
        }

        public override ddouble Entropy {
            get {
                ddouble k = Log1p(C);
                return k / 2 - Log(C / k);
            }
        }

        public override string ToString() {
            return $"{typeof(BradfordDistribution).Name}[c={C}]";
        }
    }
}
