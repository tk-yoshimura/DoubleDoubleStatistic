using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic.InternalUtils {
    public class QuantileCubicApprox {

        private readonly bool logscale;

        private readonly ReadOnlyCollection<double> value_table, grad_table;

        public int Samples { get; }

        public QuantileCubicApprox(ContinuousDistribution dist, int samples, bool logscale) {
            if (!int.IsPow2(samples) || samples < 16) {
                throw new ArgumentException(null, nameof(samples));
            }

            this.Samples = samples;

            List<ddouble> value_table = [], grad_table = [0];

            for (int i = 0; i <= samples; i++) {
                ddouble y = dist.Quantile(i / (ddouble)samples);

                if (logscale) {
                    y = ddouble.Log2(y + 1d);
                }

                value_table.Add(y);
            }

            if (logscale && !ddouble.IsFinite(value_table[^1])) {
                value_table[^1] = 4d * (value_table[^2] - value_table[^3]) + value_table[^4];
            }

            for (int i = 1; i < samples; i++) {
                ddouble ym = value_table[i - 1], y0 = value_table[i], yp = value_table[i + 1];
                ddouble g = Grad(ym, y0, yp);

                grad_table.Add(g);
            }

            grad_table[0] = ddouble.Max(0d, 2 * grad_table[1] - grad_table[2]);
            grad_table.Add(ddouble.Max(0d, 2 * grad_table[^1] - grad_table[^2]));

            this.logscale = logscale;
            this.value_table = new ReadOnlyCollection<double>(value_table.Select(d => (double)d).ToArray());
            this.grad_table = new ReadOnlyCollection<double>(grad_table.Select(d => (double)d).ToArray());
        }

        public double QuantileApprox(double p) {
            Debug.Assert(p >= 0d && p <= 1d);

            int index = int.Clamp((int)double.Floor(p * Samples), 0, Samples - 1);

            double v = p * Samples - index;

            double v0 = value_table[index], v1 = value_table[index + 1];
            double g0 = grad_table[index], g1 = grad_table[index + 1];

            double a = v0;
            double b = g0;
            double c = -3 * (v0 - v1) - 2 * g0 - g1;
            double d = 2 * (v0 - v1) + g0 + g1;

            double x = a + v * (b + v * (c + v * d));

            if (logscale) {
                x = double.Exp2(x) - 1;
            }

            return x;
        }

        private static ddouble Grad(ddouble vm1, ddouble v0, ddouble vp1) {
            ddouble vv = (vm1 - v0) * (v0 - vp1);

            return (vv > 0d) ? (2d * vv / (vp1 - vm1)) : 0d;
        }
    }
}
