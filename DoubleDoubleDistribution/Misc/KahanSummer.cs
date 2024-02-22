using DoubleDouble;

namespace DoubleDoubleDistribution {
    internal struct KahanSummer {
        private ddouble acc = 0d, carry = 0d;

        public KahanSummer() { }

        private KahanSummer(ddouble acc, ddouble carry) {
            this.acc = acc;
            this.carry = carry;
        }

        public void Add(ddouble v) { 
            ddouble d = v - carry;
            ddouble acc_next = acc + d;

            carry = (acc_next - acc) - d;
            acc = acc_next;
        }

        public static implicit operator ddouble(KahanSummer kahan) {
            return kahan.acc;
        }

        public static KahanSummer operator* (KahanSummer kahan, ddouble s) {
            return new(kahan.acc * s, kahan.carry * s);
        }
    }
}
