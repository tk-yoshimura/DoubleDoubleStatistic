using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic.Utils {
    internal class CDFSegmentCache {
        private readonly Func<ddouble, ddouble> f;
        private readonly ddouble xa, xb, xh;
        private readonly ReadOnlyCollection<ddouble> lower_segment_table, upper_segment_table;

        public int Samples { get; }

        public ReadOnlyCollection<ddouble> LowerSegmentTable => lower_segment_table;
        public ReadOnlyCollection<ddouble> UpperSegmentTable => upper_segment_table;

        public CDFSegmentCache(ddouble a, ddouble b, Func<ddouble, ddouble> f, int samples) {
            if (!int.IsPow2(samples)) {
                throw new ArgumentException(null, nameof(samples));
            }

            ddouble range = b - a;

            List<ddouble> table = [];

            this.Samples = samples;
            this.f = f;
            this.xa = a;
            this.xb = b;
            this.xh = range / samples;

            ddouble x0 = a, y0 = f(a);

            for (int i = 0; i < samples; i++) {
                ddouble x1 = x0 + xh, y1 = f(x1);

                ddouble eps = (y0 + y1) * xh * 5e-29;

                ddouble y = GaussKronrodIntegral.AdaptiveIntegrate(f, x0, x1, eps, 2048).value;

                table.Add(y);

                (x0, y0) = (x1, y1);
            }

            ddouble lower_sum = 0d, upper_sum = 0d;
            ddouble[] lower_segment_table = new ddouble[samples + 1], upper_segment_table = new ddouble[samples + 1];

            for (int i = 0; i < samples; i++) {
                lower_sum += table[i];
                lower_segment_table[i + 1] = lower_sum;
            }

            for (int i = samples - 1; i >= 0; i--) {
                upper_sum += table[i];
                upper_segment_table[i] = upper_sum;
            }

            this.lower_segment_table = new(lower_segment_table);
            this.upper_segment_table = new(upper_segment_table);
        }

        public ddouble Lower(ddouble x) {
            Debug.Assert(x >= xa && x <= xb);

            int index = int.Clamp((int)ddouble.Floor((x - xa) / xh), 0, Samples);

            ddouble x0 = xa + index * xh;

            ddouble eps = lower_segment_table[int.Max(1, index)] * 1e-28;

            ddouble cdf = lower_segment_table[index] +
                (eps > 0d ? GaussKronrodIntegral.AdaptiveIntegrate(f, x0, x, eps, 2048).value : 0d);

            return cdf;
        }

        public ddouble Upper(ddouble x) {
            Debug.Assert(x >= xa && x <= xb);

            int index = int.Clamp((int)ddouble.Ceiling((x - xa) / xh), 0, Samples);

            ddouble x0 = xa + index * xh;

            ddouble eps = upper_segment_table[int.Min(Samples - 1, index)] * 1e-28;

            ddouble ccdf = upper_segment_table[index] +
                (eps > 0d ? GaussKronrodIntegral.AdaptiveIntegrate(f, x, x0, eps, 2048).value : 0d);

            return ccdf;
        }
    }
}
