using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class TukeyLambdaDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<TukeyLambdaDistribution> {

        public ddouble Lambda { get; }

        private readonly ReadOnlyCollection<ddouble> quantile_table;
        private const int quantile_samples = 256;

        private readonly ddouble lambda_m1, lambda_m2;

        public TukeyLambdaDistribution(ddouble lambda) {
            ParamAssert.ValidateShape(nameof(lambda), IsFinite(lambda));

            Lambda = lambda;

            this.lambda_m1 = lambda - 1d;
            this.lambda_m2 = lambda - 2d;

            ddouble[] quantile_table = new ddouble[quantile_samples + 1];
            quantile_table[0] = lambda > 0d ? 1d / lambda : PositiveInfinity;
            for (int i = 1; i < quantile_samples; i++) {
                quantile_table[i] = -Quantile((ddouble)i / (quantile_samples * 2), Interval.Lower);
            }
            quantile_table[^1] = 0d;

            this.quantile_table = new ReadOnlyCollection<ddouble>(quantile_table.Reverse().ToArray());
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            x = Abs(x);

            if ((Lambda > 0d && x * Lambda > 1d) || IsPositiveInfinity(x)) {
                return 0d;
            }

            if (Abs(lambda_m1) < 1e-30d) {
                return 0.5d;
            }

            if (Abs(lambda_m2) < 1e-30d) {
                return 1d;
            }

            ddouble cdf = CDF(x, Interval.Upper);

            if (cdf == 0d) {
                return 0d;
            }

            ddouble pdf = 1d / (Pow(cdf, lambda_m1) + Pow(1d - cdf, lambda_m1));

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) && interval == Interval.Upper) {
                return 1d - CDF(-x, Interval.Upper);
            }
            if (IsPositive(x) && interval == Interval.Lower) {
                return 1d - CDF(-x, Interval.Lower);
            }

            x = Abs(x);

            if ((Lambda > 0d && x * Lambda > 1d) || IsPositiveInfinity(x)) {
                return 0d;
            }

            if (Abs(lambda_m1) < 1e-30d) {
                return 0.5d - 0.5d * x;
            }

            if (Abs(lambda_m2) < 1e-30d) {
                return 0.5d - x;
            }

            int index = Indexer.BisectionSearch(x, quantile_table);

            ddouble x0 = quantile_table[index], x1 = quantile_table[index + 1];
            ddouble p0 = (ddouble)(quantile_samples - index - 1) / (quantile_samples * 2);
            ddouble p1 = (ddouble)(quantile_samples - index) / (quantile_samples * 2);

            ddouble p = (x - x0) / (x1 - x0) * (p1 - p0) + p0;

            for (int i = 0; i < 8; i++) {
                ddouble q = (Abs(Lambda) >= 1e-64d)
                    ? (Pow(p, Lambda) - Pow(1d - p, Lambda)) / Lambda
                    : Log(p / (1d - p));

                ddouble dq = Pow(p, lambda_m1) + Pow(1d - p, lambda_m1);

                ddouble dp = (x + q) / dq;

                if (!IsFinite(dp)) {
                    break;
                }

                p = Clamp(p - dp, p0, p1);

                if (Abs(dp) <= Abs(p) * 1e-29) {
                    break;
                }
            }

            return p;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return -Quantile(p, Interval.Lower);
            }

            if (p == 0d) {
                return interval == Interval.Lower ? Support.min : Support.max;
            }
            if (p == 1d) {
                return interval == Interval.Lower ? Support.max : Support.min;
            }

            if (Abs(Lambda) >= 1e-64d) {
                if (Abs(lambda_m1) < 1e-30d) {
                    return 2d * p - 1d;
                }

                if (Abs(lambda_m2) < 1e-30d) {
                    return p - 0.5d;
                }

                ddouble x = (Pow(p, Lambda) - Pow(1d - p, Lambda)) / Lambda;

                return x;
            }
            else {
                return Log(p / (1d - p));
            }
        }

        public override double Sample(Random random) {
            double p = random.NextUniformOpenInterval01();
            double lambda = (double)Lambda;

            if (double.Abs(lambda) >= 1e-64d) {
                double x = (double.Pow(p, lambda) - double.Pow(1d - p, lambda)) / lambda;

                return x;
            }
            else {
                return double.Log(p / (1d - p));
            }
        }

        public override (ddouble min, ddouble max) Support => Lambda > 0d
            ? (-1d / Lambda, 1d / Lambda)
            : (NegativeInfinity, PositiveInfinity);

        public override ddouble Median => 0d;

        public override ddouble Mode => (Lambda < 1d || Lambda > 2d) ? 0d : NaN;

        public override ddouble Mean => (Lambda > -1d) ? 0d : NaN;

        public override ddouble Variance => (Abs(Lambda) < 0.03125)
            ? VariancePade.Value(Lambda)
            : (Lambda > -0.5d)
            ? 2d * (1d / (2d * Lambda + 1d) - Square(Gamma(Lambda + 1d)) / Gamma(2d * Lambda + 2d)) / Square(Lambda)
            : NaN;

        public override ddouble Skewness => (Lambda * 3d > -1d) ? 0d : NaN;

        public override ddouble Kurtosis {
            get {
                if (Lambda < -0.25d) {
                    return NaN;
                }

                if (Abs(Lambda) < 0.03125) {
                    return KurtosisPade.Value(Lambda);
                }

                ddouble g1 = Gamma(Lambda + 1d), g2 = Gamma(2d * Lambda + 1d);
                ddouble g3 = Gamma(3d * Lambda + 1d), g4 = Gamma(4d * Lambda + 1d);

                ddouble kurtosis = Square(g2 * (2d * Lambda + 1d) / (g1 * g1 - g2)) * (3d * g2 * g2 - 4d * g1 * g3 + g4)
                    / ((8d * Lambda + 2d) * g4) - 3d;

                return kurtosis;
            }
        }

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            Abs(lambda_m2) < 1e-30 ? 0d :
            Abs(Lambda) < 1e-64 ? 2d :
            GaussKronrodIntegral.AdaptiveIntegrate(
                p => Log(Pow(p, lambda_m1) + Pow(1d - p, lambda_m1)),
                0d, 1d, 1e-28, 65536).value;

        public static (TukeyLambdaDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (TukeyLambdaDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble lambda = BisectionMinimizeSearch1D.Search(
                lambda => {
                    try {
                        TukeyLambdaDistribution dist = new(lambda);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (-8d, 8d), iter: 32
            );

            try {
                TukeyLambdaDistribution dist = new(lambda);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(TukeyLambdaDistribution).Name}[lambda={Lambda}]";
        }

        public override string Formula => "F^-1(p; lambda) := (p^lambda - (1 - p)^lambda) / lambda";

        private static class VariancePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_nz = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 1, 0xD28D3312983E9918uL, 0x73D8912200BACE5EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0x95D78526AA5840C3uL, 0x523489321A24926BuL), (+1, 2, 0x858AC4BA1908C63FuL, 0xE906D79FB7BBAB95uL)),
                ((-1, -4, 0xCA986AD435C4D0C4uL, 0x517E477646181C43uL), (+1, 2, 0xC062A603749C15AFuL, 0x1A2AFF4CBC0C80F8uL)),
                ((-1, -6, 0xFD8C901C9B884880uL, 0x37FCF60AA6AB342BuL), (+1, 1, 0xFAB35215969575BAuL, 0x09AE33ABEE8298A3uL)),
                ((+1, -5, 0xA9756B73FF62C91AuL, 0xB9E84C8A288DE931uL), (+1, 0, 0xA2EA23F21FB2DD94uL, 0x2E0435CC2FBE6598uL)),
                ((-1, -7, 0xAD87E205230E67BBuL, 0x048230F07DE99B6DuL), (+1, -3, 0xD6BC75F83CC25966uL, 0x363B8C3D9BCB6661uL)),
                ((+1, -10, 0xBB4B3778C3272E3EuL, 0xF6A230BD0C2D7599uL), (+1, -7, 0xFA459F7F3886E6A9uL, 0x94F4C8C5D8CE97DBuL)),
                ((-1, -14, 0xBDBE9D2DA5BC2CFDuL, 0x7753938459A33CBDuL), Zero),
            }));

            public static ddouble Value(ddouble x) {
                Debug.Assert(Abs(x) < 0.03125);

                Debug.WriteLine("pade approx called");

                return ApproxUtil.Pade(x, pade_nz);
            }
        }

        private static class KurtosisPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_nz = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0x9999999999999999uL, 0x9999999999999999uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 2, 0xA4244B45CE08B2C4uL, 0x11A391DB46D473CEuL), (+1, 2, 0xF8C8CC3F5A5EC73CuL, 0x9A60F0884EA78773uL)),
                ((-1, 4, 0xE05F0E590CD2D6F9uL, 0xFEC395164FC8CAF6uL), (+1, 3, 0xFBE0FA3B76C4A06FuL, 0x7A9BDC504FADBD72uL)),
                ((-1, 2, 0x9ED4602981B62740uL, 0x28D6286C139B0652uL), (-1, 1, 0xD67E5FECD1F76E86uL, 0x887C126F4CD034F7uL)),
                ((+1, 5, 0xA91FE1B532C84216uL, 0x7E678493B421683AuL), (-1, 4, 0xD38F8451E02B7AE6uL, 0xCAEB4E3E5A8A6207uL)),
                ((+1, 3, 0xF465FB7BB66F2382uL, 0x981CFFC765DDC4C8uL), (-1, 3, 0xB21D009F05928E6FuL, 0x1F6B30750D7B8DB2uL)),
                ((+1, -1, 0xD172C3713407D56AuL, 0x97F5E4EBB7C4A42FuL), (-1, 0, 0x823F470CFC12C708uL, 0xDE9872EB5EF7C115uL)),
                ((+1, -1, 0x8AACADC8D2B9FF67uL, 0x045A43F1F43BDB73uL), (-1, -1, 0x9F2D45D2919205A7uL, 0x5AC286260FE19228uL)),
                ((-1, -2, 0xA1836B0B8ED9FD02uL, 0x10FE6E3A935FB99EuL), Zero),
            }));

            public static ddouble Value(ddouble x) {
                Debug.Assert(Abs(x) < 0.03125);

                Debug.WriteLine("pade approx called");

                return ApproxUtil.Pade(x, pade_nz);
            }
        }
    }
}
