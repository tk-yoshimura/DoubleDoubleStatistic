using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class IrwinHallDistribution : ContinuousDistribution,
        IAdditionOperators<IrwinHallDistribution, IrwinHallDistribution, IrwinHallDistribution> {

        public int N { get; }

        private readonly ReadOnlyCollection<ddouble> binom, cdf_intway;
        private readonly ddouble pdf_norm, cdf_norm, quantile_norm;

        public IrwinHallDistribution(int n) {
            ValidateShape(n, n => n > 0 && n <= 64);

            N = n;

            pdf_norm = TaylorSequence[n - 1];
            cdf_norm = TaylorSequence[n];
            quantile_norm = Factorial[n];

            List<ddouble> coefs = [], intway = [0];
            for (int k = 0; k <= n / 2; k++) {
                ddouble c = Binomial(n, k);
                coefs.Add(c);
            }
            for (int k = 1; k <= n / 2; k++) {
                ddouble c = 0d;
                for (int i = 0; i < k; i++) {
                    ddouble s = coefs[i] * Pow(k - i, N);
                    if (i % 2 == 0) {
                        c += s;
                    }
                    else {
                        c -= s;
                    }
                }
                intway.Add(c);
            }

            if (n % 2 == 0) {
                intway.Add(2 * intway[^1] - intway[^2]);
            }
            else {
                intway.Add(quantile_norm - intway[^1]);
            }

            binom = Array.AsReadOnly(coefs.ToArray());
            cdf_intway = Array.AsReadOnly(intway.ToArray());
        }

        public override ddouble PDF(ddouble x) {
            if (IsNegative(x) || x > N) {
                return 0d;
            }
            if (IsNaN(x)) {
                return NaN;
            }

            x = (x * 2 < N) ? x : N - x;

            if (x <= 1d) {
                ddouble pdf = Pow(x, N - 1) * pdf_norm;
                return pdf;
            }
            else {
                ddouble pdf = 0d;
                for (int k = 0, xn = (int)Floor(x); k <= xn; k++) {
                    ddouble xk = x - k, xpownm1 = Pow(xk, N - 1);

                    if (k % 2 == 0) {
                        pdf += binom[k] * xpownm1;
                    }
                    else {
                        pdf -= binom[k] * xpownm1;
                    }
                }

                pdf *= pdf_norm;
                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                x = N - x;
            }

            if (x * 2 > N) {
                return 1d - CDF(x: N - x, interval: Interval.Lower);
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (x <= 1d) {
                ddouble cdf = Pow(x, N) * cdf_norm;
                return cdf;
            }
            else {
                ddouble cdf = 0d;
                for (int k = 0, xn = (int)Floor(x); k <= xn; k++) {
                    ddouble xk = x - k, xpown = Pow(xk, N);

                    if (k % 2 == 0) {
                        cdf += binom[k] * xpown;
                    }
                    else {
                        cdf -= binom[k] * xpown;
                    }
                }

                cdf *= cdf_norm;
                return cdf;
            }
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

            p *= quantile_norm;

            if (p <= 1d) {
                ddouble x = RootN(p, N);
                return x;
            }
            else {
                int cdf_idx;
                for (cdf_idx = cdf_intway.Count - 2; cdf_idx > 0; cdf_idx--) {
                    if (p >= cdf_intway[cdf_idx]) {
                        break;
                    }
                }

                Debug.Assert(cdf_idx >= 0 && cdf_idx < cdf_intway.Count, nameof(cdf_idx));
                Debug.Assert(p >= cdf_intway[cdf_idx], nameof(cdf_idx));

                ddouble x = cdf_idx + (p - cdf_intway[cdf_idx]) / (cdf_intway[cdf_idx + 1] - cdf_intway[cdf_idx]);

                for (int iter = 0; iter < 16; iter++) {
                    ddouble pdf = 0d, cdf = 0d;
                    for (int k = 0; k <= cdf_idx; k++) {
                        ddouble xk = x - k, xpownm1 = Pow(xk, N - 1), xpown = xpownm1 * xk;

                        if (k % 2 == 0) {
                            pdf += binom[k] * xpownm1;
                            cdf += binom[k] * xpown;
                        }
                        else {
                            pdf -= binom[k] * xpownm1;
                            cdf -= binom[k] * xpown;
                        }
                    }

                    pdf *= N;

                    ddouble dx = (cdf - p) / pdf;

                    x -= dx;
                    x = Clamp(x, cdf_idx, cdf_idx + 1);

                    if (Abs(dx) < 1e-30) {
                        break;
                    }
                }

                return x;
            }
        }

        public override bool AdditiveClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, N);

        public override ddouble Mean => N / 0.5d;

        public override ddouble Median => N / 0.5d;

        public override ddouble Mode => N / 0.5d;

        public override ddouble Variance => (ddouble)N / 12d;

        public override ddouble Skewness => 0;

        public override ddouble Kurtosis => -(ddouble)6 / (5 * N);

        public override ddouble Entropy => throw new NotImplementedException();

        public static IrwinHallDistribution operator +(IrwinHallDistribution dist1, IrwinHallDistribution dist2) {
            return new(dist1.N + dist2.N);
        }

        public override string ToString() {
            return $"{typeof(IrwinHallDistribution).Name}[n={N}]";
        }
    }
}
