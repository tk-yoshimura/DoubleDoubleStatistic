using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class IrwinHallDistribution : ContinuousDistribution,
        IAdditionOperators<IrwinHallDistribution, IrwinHallDistribution, IrwinHallDistribution> {

        public int N { get; }

        private readonly ReadOnlyCollection<ddouble> cdf_intway;
        private readonly List<ReadOnlyCollection<ddouble>> pdf_table = [], cdf_table = [];
        private readonly ddouble pdf_norm, cdf_norm, quantile_norm;

        public IrwinHallDistribution(int n) {
            ValidateShape(n, n => n > 0 && n <= 256);

            N = n;

            pdf_norm = TaylorSequence[n - 1];
            cdf_norm = TaylorSequence[n];
            quantile_norm = Factorial[n];

            for (int k = 0; k < (n + 1) / 2; k++) {
                Polynomial pdf = PolynomialCoef.PDF(n, k);
                Polynomial cdf = PolynomialCoef.CDF(n, k);

                pdf_table.Add(new ReadOnlyCollection<ddouble>(pdf.Coef.Select(c => (ddouble)c).Reverse().ToArray()));
                cdf_table.Add(new ReadOnlyCollection<ddouble>(cdf.Coef.Select(c => (ddouble)c).Reverse().ToArray()));
            }

            List<ddouble> intway = [0];

            for (int k = 1; k < cdf_table.Count; k++) {
                intway.Add(cdf_table[k][^1]);
            }

            if (n % 2 != 0) {
                intway.Add(quantile_norm - intway[^1]);
            }
            else {
                intway.Add(quantile_norm / 2d);
            }

            cdf_intway = Array.AsReadOnly(intway.ToArray());
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || x > N) {
                return 0d;
            }

            if (x * 2d > N) {
                x = N - x;
            }

            if (x <= 1d) {
                ddouble pdf = Pow(x, N - 1) * pdf_norm;
                return pdf;
            }

            int index = int.Clamp((int)Floor(x), 0, pdf_table.Count - 1);
            ddouble c = x - index;

            ReadOnlyCollection<ddouble> coef = pdf_table[index];

            ddouble y = coef[0];

            for (int i = 1; i < coef.Count; i++) {
                y = coef[i] + y * c;
            }

            y *= pdf_norm;

            return y;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                x = N - x;
            }

            if (x * 2d > N) {
                return 1d - CDF(x: N - x, interval: Interval.Lower);
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (x <= 1d) {
                ddouble cdf = Pow(x, N) * cdf_norm;
                return cdf;
            }

            int index = int.Clamp((int)Floor(x), 0, cdf_table.Count - 1);
            ddouble c = x - index;

            ReadOnlyCollection<ddouble> coef = cdf_table[index];

            ddouble y = coef[0];

            for (int i = 1; i < coef.Count; i++) {
                y = coef[i] + y * c;
            }

            y *= cdf_norm;

            return y;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                p = 1d - p;
            }

            if (p > 0.5d) {
                return N - Quantile(1d - p, interval: Interval.Lower);
            }

            ddouble v = p * quantile_norm;

            if (v <= 1d) {
                ddouble x = RootN(v, N);
                return x;
            }
            else {
                int index;
                for (index = cdf_intway.Count - 2; index > 0; index--) {
                    if (v >= cdf_intway[index]) {
                        break;
                    }
                }

                Debug.Assert(index >= 0 && index < cdf_intway.Count, nameof(index));
                Debug.Assert(v >= cdf_intway[index], nameof(index));

                ddouble c = (v - cdf_intway[index]) / (cdf_intway[index + 1] - cdf_intway[index]);

                ReadOnlyCollection<ddouble> pdf_coef = pdf_table[index];
                ReadOnlyCollection<ddouble> cdf_coef = cdf_table[index];

                for (int iter = 0; iter < 16; iter++) {
                    ddouble pdf = pdf_coef[0], cdf = cdf_coef[0];

                    for (int i = 1; i < pdf_coef.Count; i++) {
                        pdf = pdf_coef[i] + pdf * c;
                    }
                    for (int i = 1; i < cdf_coef.Count; i++) {
                        cdf = cdf_coef[i] + cdf * c;
                    }

                    pdf *= N;

                    ddouble dc = (cdf - v) / pdf;

                    c -= dc;
                    c = Clamp(c, 0d, 1d);

                    if (Abs(dc) < 1e-29) {
                        break;
                    }
                }

                ddouble x = index + c;

                return x;
            }
        }

        public override bool AdditiveClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, N);

        public override ddouble Mean => N * 0.5d;

        public override ddouble Median => N * 0.5d;

        public override ddouble Mode => (N > 1) ? (N * 0.5d) : NaN;

        public override ddouble Variance => (ddouble)N / 12d;

        public override ddouble Skewness => 0;

        public override ddouble Kurtosis => -(ddouble)6 / (5 * N);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            (N == 1) ? 0d :
            (N == 2) ? 0.5d :
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public static IrwinHallDistribution operator +(IrwinHallDistribution dist1, IrwinHallDistribution dist2) {
            return new(dist1.N + dist2.N);
        }

        public override string ToString() {
            return $"{typeof(IrwinHallDistribution).Name}[n={N}]";
        }

        public override string Formula => "p(x; n) := sum((-1)^k * binom(n, k) * (x - k)^(n - 1), k, 0, floor(x)) / (n - 1)!";

        private static class PolynomialCoef {
            static readonly Dictionary<(int n, int k, int j), BigInteger> table = [];

            private static BigInteger ACoef(int n, int k, int j) {
                if (n <= 0 || k < 0 || j < 0 || j >= n || k >= n) {
                    return BigInteger.Zero;
                }

                if (k == 0) {
                    return (j < n - 1) ? BigInteger.Zero : BigInteger.One;
                }

                if (table.TryGetValue((n, k, j), out BigInteger value)) {
                    return value;
                }

                BigInteger c = Binom.Value(n, k) * Binom.Value(n - 1, j) * BigInteger.Pow(k, n - j - 1);

                int s = n + k - j - 1;

                BigInteger new_value = ACoef(n, k - 1, j) + (((s & 1) == 0) ? c : -c);

                table.Add((n, k, j), new_value);

                return new_value;
            }

            public static Polynomial PDF(int n, int k) {
                if (n <= 0 || k < 0 || k >= n) {
                    return Polynomial.Zero;
                }

                List<BigInteger> coef = [];

                for (int j = 0; j < n; j++) {
                    coef.Add(ACoef(n, k, j));
                }

                Polynomial p = new([.. coef]);
                Polynomial p_sft = Polynomial.XShift(p, k);

                return p_sft;
            }

            public static Polynomial CDF(int n, int k) {
                if (n <= 0 || k < 0 || k >= n) {
                    return Polynomial.Zero;
                }

                List<BigInteger> coef = [.. Polynomial.XMul(n * PDF(n, k), 1).Coef];

                for (int i = 1; i < coef.Count; i++) {
                    if ((coef[i] % i) != 0) {
                        throw new ArithmeticException();
                    }

                    coef[i] /= i;
                }

                BigInteger c0 = BigInteger.Zero;

                for (int i = 0; i <= k; i++) {
                    c0 += (((i & 1) == 0) ? +1 : -1) * Binom.Value(n, i) * BigInteger.Pow(k - i, n);
                }

                coef[0] = c0;

                Polynomial p = new([.. coef]);

                return p;
            }
        }
    }
}
