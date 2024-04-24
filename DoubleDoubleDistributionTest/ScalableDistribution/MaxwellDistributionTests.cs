using DoubleDouble;
using DoubleDoubleDistribution;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleDistributionTest.ScalableDistribution {
    [TestClass()]
    public class MaxwellDistributionTests {
        readonly MaxwellDistribution dist1 = new(sigma: 1);
        readonly MaxwellDistribution dist2 = new(sigma: 2);

        MaxwellDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Sigma={dist.Sigma}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode - 1e-4), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode + 1e-4), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (MaxwellDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist1 = [
                0.000000000000000000e+00,
                1.236992771702342259e-02,
                4.833351460035615471e-02,
                1.045842451508725290e-01,
                1.760326633821497622e-01,
                2.563757566799805265e-01,
                3.387796111741550042e-01,
                4.165842162671447579e-01,
                4.839414490382867307e-01,
                5.363127596197392322e-01,
                5.707783918406934554e-01,
                5.861401287641712088e-01,
                5.828291804965128886e-01,
                5.626552285646524432e-01,
                5.284485778123830801e-01,
                4.836535019064274188e-01,
                4.319277321055045027e-01,
                3.767926480963080249e-01,
                3.213639748361326043e-01,
                2.681767562374651193e-01,
                2.191037561696067282e-01,
                1.753551276313331797e-01,
                1.375401328365646880e-01,
                1.057689263866487817e-01,
                7.977327141488413376e-02,
                5.902847724997570267e-02,
                4.286364021045760420e-02,
                3.055236433670266427e-02,
                2.138072602862112204e-02,
                1.469283563116408411e-02,
                9.916753566584403032e-03,
                6.574790746618506231e-03,
                4.282567224476331764e-03,
                2.740892139907563380e-03,
                1.723835745051028635e-03,
                1.065521622837132176e-03,
                6.473415148296716862e-04,
                3.865879693784790485e-04,
                2.269565163977341172e-04,
                1.309933965897971194e-04,
                7.433597573671488093e-05,
                4.147819035033313670e-05,
                2.275819632482286526e-05,
                1.227940972229913055e-05,
                6.515704825738682364e-06,
                3.400248786688290305e-06,
                1.745199568361041605e-06,
                8.810145176128747649e-07,
                4.374635651872766711e-07,
                2.136671642778333872e-07,
                1.026563920434284202e-07,
                4.851777212506460478e-08,
                2.255775339474610282e-08,
                1.031771920176008977e-08,
                4.642762214996306478e-09,
                2.055348194703086274e-09,
                8.952026000197301930e-10,
                3.836139527769571404e-10,
                1.617388020527540334e-10,
                6.709496928275307139e-11,
                2.738610599657636047e-11,
                1.099875306087707922e-11,
                4.346476209398911303e-12,
                1.690126930141479059e-12,
            ];
            ddouble[] expected_dist2 = [
                0.000000000000000000e+00,
                1.555327565183150943e-03,
                6.184963858511711293e-03,
                1.378093001669198313e-02,
                2.416675730017807736e-02,
                3.710259826158474050e-02,
                5.229212257543626452e-02,
                6.939095130852272753e-02,
                8.801633169107488108e-02,
                1.077577153902589263e-01,
                1.281878783399902633e-01,
                1.488742101039757515e-01,
                1.693898055870775021e-01,
                1.893240110106408081e-01,
                2.082921081335723790e-01,
                2.259438638789660336e-01,
                2.419707245191433653e-01,
                2.561114917733861573e-01,
                2.681563798098696161e-01,
                2.779494144242578568e-01,
                2.853891959203467277e-01,
                2.904281030175342448e-01,
                2.930700643820856044e-01,
                2.933670654826524005e-01,
                2.914145902482564443e-01,
                2.873462188093751091e-01,
                2.813276142823262216e-01,
                2.735501334368267146e-01,
                2.642242889061915401e-01,
                2.535732754445184089e-01,
                2.418267509532137094e-01,
                2.292150361111904289e-01,
                2.159638660527522513e-01,
                2.022897952454552895e-01,
                1.883963240481540125e-01,
                1.744707837532631922e-01,
                1.606819874180663021e-01,
                1.471786274186571153e-01,
                1.340883781187325596e-01,
                1.215176437828144745e-01,
                1.095518780848033641e-01,
                9.825639225047400060e-02,
                8.767756381566658985e-02,
                7.784435681172226229e-02,
                6.877006641828234401e-02,
                6.045420618344322233e-02,
                5.288446319332439083e-02,
                4.603865546122722086e-02,
                3.988663570744206688e-02,
                3.439209606817358184e-02,
                2.951423862498785133e-02,
                2.520928658087766977e-02,
                2.143182010522880210e-02,
                1.813592909194991845e-02,
                1.527618216835133214e-02,
                1.280841717703836910e-02,
                1.069036301431056102e-02,
                8.882106185807097862e-03,
                7.346417815582042057e-03,
                6.048958231705831717e-03,
                4.958376783292201516e-03,
                4.046324364228536669e-03,
                3.287395373309253115e-03,
                2.659014669265068621e-03,
            ];

            foreach ((MaxwellDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist1 = [
                0.000000000000000000e+00,
                5.170279240384146680e-04,
                4.078592964422843477e-03,
                1.344821294312090708e-02,
                3.085959578372675718e-02,
                5.782773121463019583e-02,
                9.503914701405698218e-02,
                1.423298471325276648e-01,
                1.987480430987991487e-01,
                2.626885127310585322e-01,
                3.320777391937346223e-01,
                4.045848252870223871e-01,
                4.778328104646085706e-01,
                5.495880697451595021e-01,
                6.179110703795758219e-01,
                6.812587421193831982e-01,
                7.385358700508888319e-01,
                7.890991997636636057e-01,
                8.327226214294070061e-01,
                8.695345208664325698e-01,
                8.999391668806050459e-01,
                9.245331498906097512e-01,
                9.440258797711400529e-01,
                9.591705332334475864e-01,
                9.707091134651117947e-01,
                9.793328366817909414e-01,
                9.856571375427699566e-01,
                9.902092906953897655e-01,
                9.934259629766086785e-01,
                9.956578425120543585e-01,
                9.971786978118425404e-01,
                9.981966551078220728e-01,
                9.988660157102147163e-01,
                9.992984678237238905e-01,
                9.995730145378498577e-01,
                9.997443089526717452e-01,
                9.998493509837883542e-01,
                9.999126674280152116e-01,
                9.999501855142732820e-01,
                9.999720414751227526e-01,
                9.999845595017089472e-01,
                9.999916091568193943e-01,
                9.999955130062515574e-01,
                9.999976388759741486e-01,
                9.999987773472703845e-01,
                9.999993769605293181e-01,
                9.999996875626866455e-01,
                9.999998458116485800e-01,
                9.999999251162304814e-01,
                9.999999642087289420e-01,
                9.999999831645245418e-01,
                9.999999922064061186e-01,
                9.999999964492564120e-01,
                9.999999984078835613e-01,
                9.999999992973989116e-01,
                9.999999996948413239e-01,
                9.999999998695542347e-01,
                9.999999999451186783e-01,
                9.999999999772744008e-01,
                9.999999999907380754e-01,
                9.999999999962847497e-01,
                9.999999999985331733e-01,
                9.999999999994300115e-01,
                9.999999999997819522e-01,
            ];
            ddouble[] expected_dist2 = [
                0.000000000000000000e+00,
                6.485597263361922232e-05,
                5.170279240384146680e-04,
                1.734789467079922735e-03,
                4.078592964422843477e-03,
                7.882808066171824793e-03,
                1.344821294312090708e-02,
                2.103544650785202363e-02,
                3.085959578372675718e-02,
                4.308605275746005964e-02,
                5.782773121463019583e-02,
                7.514368802791128454e-02,
                9.503914701405698218e-02,
                1.174668786283585981e-01,
                1.423298471325276648e-01,
                1.694849999637876303e-01,
                1.987480430987991487e-01,
                2.298990221873312700e-01,
                2.626885127310585322e-01,
                2.968442138569268307e-01,
                3.320777391937346223e-01,
                3.680914045966626191e-01,
                4.045848252870223871e-01,
                4.412611534601585395e-01,
                4.778328104646085706e-01,
                5.140265940586649096e-01,
                5.495880697451595021e-01,
                5.842851845217601481e-01,
                6.179110703795758219e-01,
                6.502860324695791672e-01,
                6.812587421193831982e-01,
                7.107066770875858364e-01,
                7.385358700508888319e-01,
                7.646800409931069353e-01,
                7.890991997636636057e-01,
                8.117778116365563790e-01,
                8.327226214294070061e-01,
                8.519602309704510379e-01,
                8.695345208664325698e-01,
                8.855040011359291974e-01,
                8.999391668806050459e-01,
                9.129199253330086972e-01,
                9.245331498906097512e-01,
                9.348704056322240108e-01,
                9.440258797711400529e-01,
                9.520945399191581560e-01,
                9.591705332334475864e-01,
                9.653458307328354193e-01,
                9.707091134651117947e-01,
                9.753448908740720569e-01,
                9.793328366817909414e-01,
                9.827473238401525091e-01,
                9.856571375427699566e-01,
                9.881253438161813341e-01,
                9.902092906953897655e-01,
                9.919607192865828038e-01,
                9.934259629766086785e-01,
                9.946462145135283039e-01,
                9.956578425120543585e-01,
                9.964927410010230302e-01,
                9.971786978118425404e-01,
                9.977397698091440281e-01,
                9.981966551078220728e-01,
                9.985670544423983497e-01,
            ];

            foreach ((MaxwellDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}