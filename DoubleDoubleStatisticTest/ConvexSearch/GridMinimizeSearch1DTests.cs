using DoubleDouble;
using DoubleDoubleStatistic.Optimizer;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ConvexSearch {
    [TestClass()]
    public class GridMinimizeSearch1DTests {
        [TestMethod()]
        public void SearchTest1() {
            for (ddouble p = 0.5; p <= 4; p += 0.125) {
                ddouble f(ddouble x) => -ddouble.Pow(x, p) * ddouble.Exp(-x);

                ddouble x = GridMinimizeSearch1D.Search(f, (0, 5), 1000);

                Assert.IsTrue(f(x) <= f(x + 1e-15), $"{p},{x}");
                Assert.IsTrue(f(x) <= f(x - 1e-15), $"{p},{x}");
            }
        }

        [TestMethod()]
        public void SearchTest2() {
            for (ddouble p = 0.5; p <= 4; p += 0.125) {
                ddouble f(ddouble x) => -ddouble.Pow(x, p) * ddouble.Exp(-x);

                ddouble x = GridMinimizeSearch1D.Search(f, (0, 1), 1000);

                if (p <= 1) {
                    Assert.IsTrue(f(x) <= f(x + 1e-15), $"{p},{x}");
                    Assert.IsTrue(f(x) <= f(x - 1e-15), $"{p},{x}");
                }
                else {
                    Assert.AreEqual(1, x);
                }

            }
        }
    }
}
