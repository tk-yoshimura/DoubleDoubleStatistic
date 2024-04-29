using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic {
    internal class QuantileBuilder {

        private readonly ddouble xa, xb, xh;
        private readonly ReadOnlyCollection<ddouble> cdf_table, cdf_log2_table, weight_table;

        public int Samples { get; }

        public QuantileBuilder(ddouble a, ddouble b, ReadOnlyCollection<ddouble> cdf_table, int samples) {
            if (!int.IsPow2(samples) || cdf_table.Count != samples + 1) {
                throw new ArgumentException(null, nameof(samples));
            }

            ddouble range = b - a;

            List<ddouble> log2_table = [], weights = [0];

            this.Samples = samples;
            this.xa = a;
            this.xb = b;
            this.xh = range / samples;

            for (int i = 0; i <= samples; i++) {
                ddouble y = cdf_table[i];

                log2_table.Add(ddouble.Log2(y));
            }

            for (int i = 1; i < samples; i++) {
                ddouble ym = cdf_table[i - 1], y0 = cdf_table[i], yp = cdf_table[i + 1];
                ddouble log2_ym = log2_table[i - 1], log2_yp = log2_table[i + 1];

                ddouble yc = (ym + yp) / 2;
                ddouble log2_yc = (log2_ym + log2_yp) / 2;

                ddouble err_linear = ddouble.Abs(y0 - yc);
                ddouble err_log2 = ddouble.Abs(y0 - ddouble.Pow2(log2_yc));

                ddouble weight = err_log2 / (err_linear + err_log2);

                weights.Add(weight);
            }

            weights[0] = weights[1];
            weights.Add(weights[^1]);

            this.cdf_table = cdf_table;
            cdf_log2_table = new(log2_table);
            weight_table = new(weights);
        }

        public QuantileBuilder(ddouble a, ddouble b, Func<ddouble, ddouble> f, int samples) {
            if (!int.IsPow2(samples)) {
                throw new ArgumentException(null, nameof(samples));
            }

            ddouble range = b - a;

            List<ddouble> table = [], log2_table = [], weights = [0];

            this.Samples = samples;
            this.xa = a;
            this.xb = b;
            this.xh = range / samples;

            for (int i = 0; i <= samples; i++) {
                ddouble x = a + xh * i;

                ddouble y = f(x);

                table.Add(y);
                log2_table.Add(ddouble.Log2(y));
            }

            for (int i = 1; i < samples; i++) {
                ddouble ym = table[i - 1], y0 = table[i], yp = table[i + 1];
                ddouble log2_ym = log2_table[i - 1], log2_yp = log2_table[i + 1];

                ddouble yc = (ym + yp) / 2;
                ddouble log2_yc = (log2_ym + log2_yp) / 2;

                ddouble err_linear = ddouble.Abs(y0 - yc);
                ddouble err_log2 = ddouble.Abs(y0 - ddouble.Pow2(log2_yc));

                ddouble weight = err_log2 / (err_linear + err_log2);

                weights.Add(weight);
            }

            weights[0] = weights[1];
            weights.Add(weights[^1]);

            cdf_table = new(table);
            cdf_log2_table = new(log2_table);
            weight_table = new(weights);
        }

        public (ddouble x, ddouble x0, ddouble x1) Estimate(ddouble p) {
            if (p < cdf_table[0]) {
                if (xa < xb) {
                    return (xa, ddouble.NegativeInfinity, xa);
                }
                else {
                    return (xa, ddouble.PositiveInfinity, xa);
                }
            }
            if (p > cdf_table[^1]) {
                if (xa < xb) {
                    return (xb, xb, ddouble.PositiveInfinity);
                }
                else {
                    return (xb, xb, ddouble.NegativeInfinity);
                }
            }

            int index = 0;

            for (int h = Samples / 2; h >= 1; h /= 2) {
                for (int i = index; i <= Samples - h; i += h) {
                    if (cdf_table[i + h] > p) {
                        index = i;
                        break;
                    }
                }
            }

            Debug.Assert(index >= 0 && index < Samples);
            Debug.Assert(cdf_table[index] <= p && cdf_table[index + 1] >= p);

            ddouble weight = (weight_table[index] + weight_table[index + 1]) / 2;

            ddouble c =
                weight * (p - cdf_table[index]) / (cdf_table[index + 1] - cdf_table[index]) +
                (1d - weight) * (ddouble.Log2(p) - cdf_log2_table[index]) / (cdf_log2_table[index + 1] - cdf_log2_table[index]);

            if (ddouble.IsNaN(c)) {
                c = 0d;
            }

            ddouble x0 = xa + xh * index, x1 = x0 + xh;
            ddouble x = xa + xh * (index + c);

            return (x, x0, x1);
        }
    }
}
