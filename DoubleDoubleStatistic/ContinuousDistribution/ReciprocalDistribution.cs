using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class ReciprocalDistribution : ContinuousDistribution {

        public ddouble A { get; }
        public ddouble B { get; }

        private readonly ddouble pdf_norm;

        public ReciprocalDistribution(ddouble a, ddouble b) {
            ValidateShape(a, a => a > 0d);
            ValidateShape(b, b => b > a);

            A = a;
            B = b;
            pdf_norm = 1d / Log(b / a);
        }

        public override ddouble PDF(ddouble x) {
            if (x < A || x > B) {
                return 0d;
            }

            ddouble pdf = pdf_norm / x;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                if (x < A) {
                    return 0d;
                }

                if (x > B) {
                    return 1d;
                }

                ddouble cdf = Log(x / A) * pdf_norm;

                return cdf;
            }
            else {
                if (x < A) {
                    return 1d;
                }

                if (x > B) {
                    return 0d;
                }

                ddouble cdf = 1d - Log(x / A) * pdf_norm;

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

            ddouble x = A * Exp(p / pdf_norm);

            x = Clamp(x, A, B);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mean => (B - A) * pdf_norm;

        public override ddouble Median => Sqrt(A * B);

        public override ddouble Mode => A;

        public override ddouble Variance =>
            (B * B - A * A) * pdf_norm / 2 - Square(Mean);

        public override ddouble Skewness {
            get {
                ddouble mu = Mean;
                ddouble logba = Log(B / A);
                return (-2d * (A * A * A - B * B * B) +
                        mu * (9d * (A * A - B * B) +
                        mu * (-18d * (A - B) +
                        mu * (-6d * logba)))) /
                        (6d * logba * Cube(Sqrt(Variance)));
            }
        }

        public override ddouble Kurtosis {
            get {
                ddouble mu = Mean;
                ddouble logba = Log(B / A);
                return (-3d * (A * A * A * A - B * B * B * B) +
                        mu * (16d * (A * A * A - B * B * B) +
                        mu * (-36d * (A * A - B * B) +
                        mu * (48d * (A - B) +
                        mu * (12d * logba))))) /
                        (12d * logba * Square(Variance)) - 3d;
            }
        }

        public override ddouble Entropy {
            get {
                ddouble logba = Log(B / A);

                return -(Square(Log(A * logba)) - Square(Log(B * logba))) / (2 * logba);
            }
        }

        public override string ToString() {
            return $"{typeof(ReciprocalDistribution).Name}[a={A},b={B}]";
        }

        public override string Formula => "p(x; a, b) := 1 / (x * log(b / a))";
    }
}
