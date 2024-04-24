using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class InverseChiDistributionTests {
        readonly InverseChiDistribution dist_nu1 = new(nu: 1);
        readonly InverseChiDistribution dist_nu2 = new(nu: 2);
        readonly InverseChiDistribution dist_nu3 = new(nu: 3);
        readonly InverseChiDistribution dist_nu4 = new(nu: 4);
        readonly InverseChiDistribution dist_nu5 = new(nu: 5);
        readonly InverseChiDistribution dist_nu8 = new(nu: 8);
        readonly InverseChiDistribution dist_nu16 = new(nu: 16);
        readonly InverseChiDistribution dist_nu32 = new(nu: 32);
        readonly InverseChiDistribution dist_nu64 = new(nu: 64);
        readonly InverseChiDistribution dist_nu128 = new(nu: 128);


        InverseChiDistribution[] Dists => [
            dist_nu1,
            dist_nu2,
            dist_nu3,
            dist_nu4,
            dist_nu5,
            dist_nu8,
            dist_nu16,
            dist_nu32,
            dist_nu64,
            dist_nu128,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.Nu}");
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
            foreach (InverseChiDistribution dist in Dists) {
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
            foreach (InverseChiDistribution dist in Dists) {
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
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);
                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);

                Console.WriteLine($"cdf(median)={cdf_median}");

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (InverseChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);
                ddouble ccdf_median = dist.CDF(dist.Median, Interval.Upper);

                Console.WriteLine($"ccdf(median)={ccdf_median}");

                Assert.IsTrue(ddouble.Abs(cdf_median + ccdf_median - 1) < 1e-28);

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (InverseChiDistribution dist in Dists) {
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
            foreach (InverseChiDistribution dist in Dists) {
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
            ddouble[] expected_nu1 = [
                0,
                0.0,
                3.576558093841976e-220,
                5.419202701601473e-97,
                5.253954932694547e-54,
                3.459043656483763e-34,
                1.784296331728179e-23,
                4.702915941929588e-17,
                6.466906986927222e-13,
                4.218273714909452e-10,
                4.168204518723632e-8,
                1.204512585475054e-6,
                1.511140641292466e-5,
                1.055561946532257e-4,
                4.831809015838607e-4,
                0.00161835057766494,
                0.004282567224476332,
                0.009457406722157535,
                0.01813781386731284,
                0.031115993293679,
                0.04882612636599998,
                0.07128559603767445,
                0.09812480880949784,
                0.128674516067114,
                0.1620762456717858,
                0.1973894031262017,
                0.2336796635285548,
                0.2700827827230149,
                0.305844307896537,
                0.3403390198680107,
                0.3730751128027336,
                0.4036879831240426,
                0.4319277321055045,
                0.4576435216504054,
                0.4807670004752047,
                0.5012962458358016,
                0.519281072750022,
                0.5348101364733588,
                0.5479999661223736,
                0.558985885646511,
                0.5679146735588125,
                0.5749387610026322,
                0.5802117505477805,
                0.5838850422705103,
                0.5861053699129043,
                0.5870130720465893,
                0.5867409472258631,
                0.5854135657148252,
                0.5831469321813106,
                0.5800484131303024,
                0.5762168595993801,
                0.571742869815992,
                0.5667091483199781,
                0.5611909277503498,
                0.5552564273682488,
                0.5489673287161057,
                0.5423792538519523,
                0.5355422355742728,
                0.5285011721625787,
                0.5212962615681013,
                0.513963411836301,
                0.5065346259428453,
                0.4990383602709275,
                0.4914998567262954,
                0.4839414490382867,
            ];
            ddouble[] expected_nu2 = [
                0,
                0.0,
                1.434416263021544e-218,
                1.448952183231035e-95,
                1.053576959034333e-52,
                5.549143444816128e-33,
                2.385369405562788e-22,
                5.389011233487354e-16,
                6.484052761136218e-12,
                3.759517924951032e-9,
                3.343404576344754e-7,
                8.783317247711009e-6,
                1.010098095575391e-4,
                6.512988112744107e-4,
                0.002768354079249838,
                0.00865408707464608,
                0.02146960818576076,
                0.04462344111895598,
                0.08082623480853442,
                0.1313620691926979,
                0.1958223182235548,
                0.2722841761553652,
                0.3577635202976703,
                0.448750163754141,
                0.5416865333959181,
                0.6333207795055028,
                0.7209206791469956,
                0.8023669803431469,
                0.8761577026381817,
                0.9413554870746014,
                0.9975046680666473,
                1.044537509809019,
                1.082682265892902,
                1.112380306492713,
                1.134215677082996,
                1.148857914311674,
                1.157017439498598,
                1.15941207320509,
                1.156742913329829,
                1.14967780601534,
                1.138840782656283,
                1.124806052300888,
                1.108095374625688,
                1.089177864996747,
                1.068471485190362,
                1.046345645489767,
                1.023124485704434,
                0.9990905165435465,
                0.974488392313333,
                0.9495286550251492,
                0.9243913424581681,
                0.8992293918741562,
                0.8741717997974835,
                0.8493265189225236,
                0.824783087710724,
                0.800614998114002,
                0.7768818133020392,
                0.7536310511997576,
                0.7308998517618461,
                0.7087164467558186,
                0.6871014507918344,
                0.6660689917209041,
                0.6456276975383735,
                0.6257815557319202,
                0.6065306597126334,
            ];
            ddouble[] expected_nu4 = [
                0,
                0.0,
                7.344211266670304e-216,
                3.297171190285733e-93,
                1.348578507563946e-50,
                4.545858309993373e-31,
                1.357010150720163e-20,
                2.252386735955531e-14,
                2.07489688356359e-10,
                9.505546555925574e-8,
                6.847292572354057e-6,
                1.486630886224144e-4,
                0.001436583958151667,
                0.007892662517692267,
                0.02892647527705953,
                0.07877142368388965,
                0.1717568654860861,
                0.3162242470990376,
                0.5109016323699953,
                0.7452341210710396,
                1.0026102693046,
                1.264485244367773,
                1.513842333821547,
                1.737316323947979,
                1.925996563185487,
                2.075265530283632,
                2.184091051616934,
                2.254111900881707,
                2.288738488524229,
                2.292385300272037,
                2.269877289111659,
                2.226027908521198,
                2.165364531785803,
                2.091969575479409,
                2.009406320645308,
                1.920702864090047,
                1.828373237726179,
                1.734460135810098,
                1.640588287049508,
                1.548021135252739,
                1.457716201800043,
                1.370376439686032,
                1.286496217252499,
                1.206401442678928,
                1.130283885160052,
                1.058230065166935,
                0.9902452489237624,
                0.9262731452608343,
                0.8662119042785181,
                0.8099269827119971,
                0.7572613877417312,
                0.7080437503107542,
                0.6620946175980941,
                0.6192312961030004,
                0.5792715238791367,
                0.5420362036818102,
                0.5073513882788827,
                0.4750496746251472,
                0.4449711344852142,
                0.4169638847905534,
                0.3908843808949102,
                0.366597499340073,
                0.3439764632046277,
                0.3229026520884284,
                0.3032653298563167,
            ];

            foreach ((InverseChiDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 64, i++) {
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
            ddouble[] expected_nu1 = [
                0,
                0.0,
                1.090416120702479e-224,
                5.569422482563291e-101,
                1.277750880107617e-57,
                1.639512342632598e-37,
                1.457619656286957e-26,
                6.082371125692384e-20,
                1.244192114854357e-15,
                1.151124166974652e-12,
                1.553769516341965e-10,
                5.949117219286293e-9,
                9.642606730228258e-8,
                8.519394455616588e-7,
                4.844104126359067e-6,
                1.984152777110947e-5,
                6.334248366623986e-5,
                1.667452332684336e-4,
                3.771812298382575e-4,
                7.560004469193409e-4,
                0.001374275875831695,
                0.002306621400566612,
                0.00362481438520609,
                0.00539237991292133,
                0.007660761135179466,
                0.01046721632711156,
                0.01383425638570892,
                0.01777027413546662,
                0.02227097895923279,
                0.0273212873976749,
                0.03289739164549068,
                0.03896880624555023,
                0.04550026389635835,
                0.0524533879436808,
                0.05978811082646412,
                0.06746383613161144,
                0.0754403596268003,
                0.08367857442187843,
                0.09214098965388828,
                0.1007920927078381,
                0.1095985833991157,
                0.1185295057610687,
                0.1275562997800073,
                0.1366527920314357,
                0.1457951409519765,
                0.1549617495684205,
                0.1641331559570925,
                0.1732919095331348,
                0.182422439451736,
                0.1915109199079658,
                0.2005451359088842,
                0.2095143521200002,
                0.2184091866192416,
                0.227221490790263,
                0.2359442361224782,
                0.2445714083313498,
                0.2530979089471155,
                0.2615194643249235,
                0.2698325418892723,
                0.2780342733285036,
                0.286122384391018,
                0.2940951308960607,
                0.3019512405520869,
                0.3096898601699836,
                0.317310507862914,
            ];
            ddouble[] expected_nu2 = [
                0,
                0.0,
                4.377491037053051e-223,
                1.492374761476057e-99,
                2.572209372642415e-56,
                2.646037790687623e-36,
                1.965483824163674e-25,
                7.051204120964674e-19,
                1.266416554909418e-14,
                1.04548971835682e-11,
                1.275407629526044e-9,
                4.459608175927487e-8,
                6.658361469857313e-7,
                5.458463624457853e-6,
                2.897782742867109e-5,
                1.114179377629491e-4,
                3.354626279025119e-4,
                8.363150261590222e-4,
                0.001798166661847583,
                0.00343708966290556,
                0.005976022895005942,
                0.009619231244563436,
                0.01453195939685666,
                0.02082803055723814,
                0.02856550078455037,
                0.03774886009129898,
                0.0483356546657089,
                0.06024547299993194,
                0.07336965136838292,
                0.08758056249337179,
                0.1027398149024943,
                0.1187050512493915,
                0.1353352832366128,
                0.1524948542573116,
                0.1700562018282703,
                0.1879016230625651,
                0.2059242464341987,
                0.2240283956301018,
                0.2421295056924224,
                0.2601537238121944,
                0.2780373004531941,
                0.2957258527016812,
                0.3131735615359038,
                0.3303423481456617,
                0.3472010612276297,
                0.3637246969042016,
                0.3798936650868485,
                0.3956931102718377,
                0.4111122905071876,
                0.4261440152551796,
                0.4407841408053244,
                0.4550311205348918,
                0.4688856064831717,
                0.4823500982575554,
                0.4954286351138359,
                0.5081265270661051,
                0.5204501210207022,
                0.5324065981477231,
                0.544003798969106,
                0.5552500729303865,
                0.5661541495171976,
                0.5767250282661534,
                0.5869718852955836,
                0.5969039942401866,
                0.6065306597126334,
            ];
            ddouble[] expected_nu4 = [
                0,
                0.0,
                2.245652902008216e-220,
                3.410905427062501e-97,
                3.318150090708715e-54,
                2.194094536038176e-34,
                1.137796747099193e-23,
                3.017627559523045e-17,
                4.179174631201078e-13,
                2.747960012816877e-10,
                2.739575588221942e-8,
                7.994124077344397e-7,
                1.013550579300502e-5,
                7.160599914451519e-5,
                3.317665548466222e-4,
                0.001125568766823037,
                0.003019163651122607,
                0.006762865799770363,
                0.01316435593179774,
                0.02293614680869666,
                0.03657326011743636,
                0.05429085389505305,
                0.07602256444801875,
                0.1014628256068104,
                0.130131725796285,
                0.1614443248384675,
                0.1947726676174424,
                0.2294947579160644,
                0.2650291488204851,
                0.3008564150337112,
                0.3365299714806147,
                0.3716789794062637,
                0.4060058497098382,
                0.4392804020249641,
                0.4713322410534411,
                0.5020424590071633,
                0.5313354012931795,
                0.5591709480409481,
                0.5855375580872154,
                0.6104461803324934,
                0.6339250450332825,
                0.6560152913293096,
                0.6767673563349578,
                0.6962380371680061,
                0.7144881342618166,
                0.7315805878966979,
                0.7475790271368796,
                0.7625466593151713,
                0.7765454376246874,
                0.7896354535069947,
                0.8018745089530463,
                0.8133178313597509,
                0.8240179001508992,
                0.8340243599989128,
                0.8433840002417979,
                0.8521407840682156,
                0.8603359143403443,
                0.8680079256351152,
                0.8751927942986926,
                0.8819240601068966,
                0.8882329545758699,
                0.8941485321331466,
                0.8996978012907334,
                0.9049058537019912,
                0.9097959895689501,
            ];

            foreach ((InverseChiDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 64, i++) {
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