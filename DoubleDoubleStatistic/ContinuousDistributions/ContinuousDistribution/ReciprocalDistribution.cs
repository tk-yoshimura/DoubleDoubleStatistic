using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class ReciprocalDistribution : ContinuousDistribution, 
        IFittableContinuousDistribution<ReciprocalDistribution> {

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

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01();

            double v = (double)A * double.Exp(u / (double)pdf_norm);

            return v;
        }

        public override (ddouble min, ddouble max) Support => (A, B);

        public override ddouble Mean => (B - A) * pdf_norm;

        public override ddouble Median => Sqrt(A * B);

        public override ddouble Mode => A;

        public override ddouble Variance =>
            (B * B - A * A) * pdf_norm / 2d - Square(Mean);

        public override ddouble Skewness {
            get {
                ddouble mu = Mean;
                ddouble logba = Log(B / A);
                return (-2d * (A * A * A - B * B * B) +
                        mu * (9d * (A * A - B * B) +
                        mu * (-18d * (A - B) +
                        mu * (-6d * logba)))) /
                        (6d * logba * ExMath.Pow3d2(Variance));
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

                return -(Square(Log(A * logba)) - Square(Log(B * logba))) / (2d * logba);
            }
        }

        public static (ReciprocalDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (ReciprocalDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            (ddouble u, ddouble v) = GridMinimizeSearch2D.Search(
                ((ddouble u, ddouble v) t) => {
                    ddouble a = t.u / (1d - t.u);
                    ddouble b = a + t.v / (1d - t.v);

                    try {
                        ReciprocalDistribution dist = new(a, b);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, ((1e-2d, 1e-2d), (100d / 101d, 100d / 101d)), iter: 64
            );

            try {
                ddouble a = u / (1d - u);
                ddouble b = a + v / (1d - v);
                ReciprocalDistribution dist = new(a, b);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(ReciprocalDistribution).Name}[a={A},b={B}]";
        }

        public override string Formula => "p(x; a, b) := 1 / (x * log(b / a))";
    }
}
