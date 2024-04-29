using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleStatistic.Utils {
    internal class CDFSegmentCache {
        private const int samples = 128;

        private readonly Func<ddouble, ddouble> f;
        private readonly ddouble xa, xb, xh;
        private readonly ReadOnlyCollection<ddouble> cdfsegment_table, ccdfsegment_table;

        public CDFSegmentCache(ddouble a, ddouble b, Func<ddouble, ddouble> f) {
            ddouble range = b - a;

            List<ddouble> table = [];

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

            ddouble cdf = 0d, ccdf = 0d;
            ddouble[] cdfsegment_table = new ddouble[samples + 1], ccdfsegment_table = new ddouble[samples + 1];

            for (int i = 0; i < samples; i++) {
                cdf += table[i];
                cdfsegment_table[i + 1] = cdf;
            }

            for (int i = samples - 1; i >= 0; i--) {
                ccdf += table[i];
                ccdfsegment_table[i] = ccdf;
            }

            this.cdfsegment_table = new(cdfsegment_table);
            this.ccdfsegment_table = new(ccdfsegment_table);
        }

        public ddouble Lower(ddouble x) {
            Debug.Assert(x >= xa && x <= xb);

            int index = int.Clamp((int)ddouble.Floor((x - xa) / xh), 0, samples);

            ddouble x0 = xa + index * xh;

            ddouble eps = cdfsegment_table[int.Max(1, index)] * 1e-28;

            ddouble cdf = cdfsegment_table[index] +
                (eps > 0d ? GaussKronrodIntegral.AdaptiveIntegrate(f, x0, x, eps, 2048).value : 0d);

            return cdf;
        }

        public ddouble Upper(ddouble x) {
            Debug.Assert(x >= xa && x <= xb);

            int index = int.Clamp((int)ddouble.Ceiling((x - xa) / xh), 0, samples);

            ddouble x0 = xa + index * xh;

            ddouble eps = ccdfsegment_table[int.Min(samples - 1, index)] * 1e-28;

            ddouble ccdf = ccdfsegment_table[index] +
                (eps > 0d ? GaussKronrodIntegral.AdaptiveIntegrate(f, x, x0, eps, 2048).value : 0d);

            return ccdf;
        }
    }
}
