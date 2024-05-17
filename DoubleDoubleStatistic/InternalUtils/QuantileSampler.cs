using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic.InternalUtils {
    internal class QuantileSampler {

        private readonly ReadOnlyCollection<double> linear_table, log2_table, weight_table;

        public int Samples { get; }

        public QuantileSampler(ContinuousDistribution dist, int samples) {
            if (!int.IsPow2(samples)) {
                throw new ArgumentException(null, nameof(samples));
            }

            this.Samples = samples;

            List<ddouble> linear_table = [], log2_table = [], weights = [0];

            for (int i = 0; i <= samples; i++) {
                ddouble y = dist.Quantile(i / (ddouble)samples);

                linear_table.Add(y);
                log2_table.Add(ddouble.Log2(y));
            }

            for (int i = 1; i < samples; i++) {
                ddouble ym = linear_table[i - 1], y0 = linear_table[i], yp = linear_table[i + 1];
                ddouble log2_ym = log2_table[i - 1], log2_yp = log2_table[i + 1];

                ddouble yc = (ym + yp) * 0.5d;
                ddouble log2_yc = (log2_ym + log2_yp) * 0.5d;

                ddouble err_linear = ddouble.Abs(y0 - yc);
                ddouble err_log2 = ddouble.Abs(y0 - ddouble.Pow2(log2_yc));

                ddouble weight = err_log2 / (err_linear + err_log2);

                weights.Add(weight);
            }

            weights[0] = weights[1] = linear_table[0] > 0 ? weights[1] : 1d;
            weights.Add(weights[^1]);

            this.linear_table = new ReadOnlyCollection<double>(linear_table.Select(d => (double)d).ToArray());
            this.log2_table = new ReadOnlyCollection<double>(log2_table.Select(d => (double)d).ToArray());

            weight_table = new(weights.Select(d => (double)d).ToArray());
        }

        public double QuantileApprox(double p) {
            Debug.Assert(p >= 0d && p <= 1d);

            int index = int.Clamp((int)double.Floor(p * Samples), 0, Samples - 1);

            double c = p * Samples - index;

            double w = (1 - c) * weight_table[index] + c * weight_table[index + 1];

            double x = w < 1
                ? (w * ((1 - c) * linear_table[index] + c * linear_table[index + 1])
                    + (1 - w) * double.Exp2((1 - c) * log2_table[index] + c * log2_table[index + 1]))
                : ((1 - c) * linear_table[index] + c * linear_table[index + 1]);

            return x;
        }
    }
}
