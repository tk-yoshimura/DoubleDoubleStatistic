using DoubleDouble;
using DoubleDoubleStatistic.Optimizer;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ConvexSearch {
    [TestClass()]
    public class GridMinimizeSearch2DTests {
        [TestMethod()]
        public void SearchTest1() {
            for (int xc = -9; xc <= 9; xc++) {
                for (int yc = -9; yc <= 9; yc++) {
                    ddouble xm = xc / (ddouble)10, ym = yc / (ddouble)10;

                    ddouble f((ddouble x, ddouble y) v) => ddouble.Square(v.x - xm) + ddouble.Pow(v.y - ym, 4);

                    (ddouble x, ddouble y) = GridMinimizeSearch2D.Search(f, ((-1, -1), (1, 1)), 1000);

                    Assert.AreEqual((double)xm, (double)x, 1e-10);
                    Assert.AreEqual((double)ym, (double)y, 1e-10);
                }
            }
        }

        [TestMethod()]
        public void SearchTest2() {
            ddouble f((ddouble x, ddouble y) v) => ddouble.Square(v.x - 4) + ddouble.Square(v.y - 5) + v.x * v.y / 4;

            (ddouble x, ddouble y) = GridMinimizeSearch2D.Search(f, ((2, 3), (5, 7)), 1000);

            Assert.AreEqual(24.0 / 7.0, (double)x, 1e-10);
            Assert.AreEqual(32.0 / 7.0, (double)y, 1e-10);
        }
    }
}
