using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LandauDistribution : StableDistribution<LandauDistribution>,
        IAdditionOperators<LandauDistribution, LandauDistribution, LandauDistribution>,
        ISubtractionOperators<LandauDistribution, LandauDistribution, LandauDistribution>,
        IAdditionOperators<LandauDistribution, ddouble, LandauDistribution>,
        ISubtractionOperators<LandauDistribution, ddouble, LandauDistribution>,
        IMultiplyOperators<LandauDistribution, ddouble, LandauDistribution>,
        IDivisionOperators<LandauDistribution, ddouble, LandauDistribution>,
        IFittableContinuousDistribution<LandauDistribution> {

        public override ddouble Mu { get; }

        public override ddouble C { get; }

        private readonly ddouble c_inv, bias;

        private static readonly ddouble mode_base = "-0.42931452986133525016556463510885028346";
        private static readonly ddouble median_base = "0.57563014394507821439627930892257517269";
        private static readonly ddouble entropy_base = "2.3726364400044818244844049010588577710";

        public LandauDistribution() : this(mu: 0d, c: 1d) { }

        public LandauDistribution(ddouble c) : this(mu: 0d, c: c) { }

        public LandauDistribution(ddouble mu, ddouble c) {
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(c), ParamAssert.IsFinitePositive(c));

            Mu = mu;
            C = c;

            c_inv = 1d / c;
            bias = -2d * RcpPi * Log(c);
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * c_inv + bias;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsInfinity(u)) {
                return 0d;
            }

            ddouble pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * c_inv + bias;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble cdf = CDFPade.Value(u, interval != Interval.Lower);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble x = Mu + C * (QuantilePade.Value(p, interval != Interval.Lower) - bias);

            return x;
        }

        public override double Sample(Random random) {
            double z = random.NextUniformOpenInterval01(), u = z - 0.5d;
            double w = random.NextUniformOpenInterval01();

            double r = 2d / double.Pi * (z * double.TanPi(u) * double.Pi - double.Log(double.Log(w) * double.CosPi(u) / (-2d * z * (double)C)));
            double v = r * (double)C + (double)Mu;

            return v;
        }

        public override ddouble Median => Mu + (median_base - bias) * C;

        public override ddouble Mode => Mu + (mode_base - bias) * C;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => entropy_base + Log(C);

        public override ddouble Alpha => 1d;

        public override ddouble Beta => 1d;

        public static LandauDistribution operator +(LandauDistribution dist1, LandauDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, dist1.C + dist2.C);
        }

        public static LandauDistribution operator -(LandauDistribution dist1, LandauDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, dist1.C + dist2.C);
        }

        public static LandauDistribution operator +(LandauDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.C);
        }

        public static LandauDistribution operator -(LandauDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.C);
        }

        public static LandauDistribution operator *(LandauDistribution dist, ddouble k) {
            return new((dist.Mu - 2d * RcpPi * dist.C * Log(k)) * k, dist.C * k);
        }

        public static LandauDistribution operator /(LandauDistribution dist, ddouble k) {
            return new((dist.Mu + 2d * RcpPi * dist.C * Log(k)) / k, dist.C / k);
        }

        public static (LandauDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (LandauDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileLinearFitter<LandauDistribution>.Fit(new LandauDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(LandauDistribution).Name}[mu={Mu},c={C}]";
        }

        public override string Formula => "p(x; mu, c) := stable_distribution(x; alpha = 1, beta = 1, mu, c)";

        private static class PDFPade {
            private static readonly ddouble pi_half = Ldexp(Pi, -1);
            private static readonly ddouble lambda_bias = (+1, 0, 0xB9CD764B54262B88uL, 0xCAD8DC6940262974uL);

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8644567CF64B3013uL, 0xF4536A34780ED6F4uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xA63C576140787BF0uL, 0xA959C4411D9EA82CuL), (+1, 1, 0xB36C485858CBD62AuL, 0x87D6D3F68C8EBD3CuL)),
                ((+1, -1, 0xB5FF59AD5E41100DuL, 0xE2D07FB9F81850A2uL), (+1, 1, 0xF51AB6BFA29ACEDDuL, 0x0E0D683436AF274BuL)),
                ((+1, -2, 0xE971BA5D1D2E95A4uL, 0xEB4789795D4492DAuL), (+1, 1, 0xD4B28ABA87B411BBuL, 0x476C43C022430A31uL)),
                ((+1, -3, 0xC3B2D43F97E3C7BAuL, 0xDC5A99D6EFF7E408uL), (+1, 1, 0x81AA080387E1437DuL, 0xCF00454C7394EE6EuL)),
                ((+1, -5, 0xE13F5A89426DAC76uL, 0xF33670EAB426FB6AuL), (+1, -1, 0xE8ADEE4A256083BBuL, 0xA374E9915007B384uL)),
                ((+1, -7, 0xB3EDC7EA72857D5DuL, 0xE9FBF369B6631E05uL), (+1, -2, 0x9CBF8F31656B4380uL, 0x979BD394C55ACE7FuL)),
                ((+1, -10, 0xC2AA21ECB9A4E171uL, 0x044F3C5B412EE6F9uL), (+1, -4, 0x9ED591DA39D156F0uL, 0xD97B65EAC4C5B882uL)),
                ((+1, -13, 0x84F97BC74A90EF82uL, 0x9BC2FFE8A4AB70A1uL), (+1, -7, 0xEE5C78F03B3C12D4uL, 0x2E0E375AE818F9ACuL)),
                ((+1, -18, 0xC0EB25DF33A2E7EBuL, 0xD8552E50ECB250F7uL), (+1, -10, 0xFECE9C1BE1B40B8DuL, 0x36BA860F53B58FA1uL)),
                ((+1, -24, 0xBB67BB8D8E13556CuL, 0xE47EC5F657865A09uL), (+1, -13, 0xB35C5B12519F0692uL, 0xD53B6A3EC2211D1FuL)),
                (Zero, (+1, -17, 0x8B7B40609137DAB1uL, 0x0D6147C506FCAC24uL)),
                (Zero, (+1, -23, 0x91DBFA4FA3088705uL, 0xD448CB4AE53CF216uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA774BBCD9C1D5B20uL, 0x7CCC9D5E31A8D9EDuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8A57AF0E51B3444AuL, 0x4C24069129A46582uL), (+1, 1, 0x8C65FB1DD75A11FAuL, 0x4E3A7C50A89A3B6EuL)),
                ((+1, -3, 0xC7EC4088C0F2610CuL, 0x68E90968B86A9B1CuL), (+1, 1, 0x903B51668999B69AuL, 0xF11D44E537E9C18DuL)),
                ((+1, -4, 0xA65860D3C57BF4F8uL, 0x733A7725900F0AF9uL), (+1, 0, 0xB5C5EE2558B19712uL, 0x575F47A33D4AAFD5uL)),
                ((+1, -6, 0xAF830A0EE73617F2uL, 0x3B87909748FA9B82uL), (+1, -1, 0x9B24837E19BB4425uL, 0x0E5A25E806383507uL)),
                ((+1, -9, 0xF121864CA7D470E6uL, 0x23D734ED86B34F4DuL), (+1, -3, 0xBAFDE4505975E6ABuL, 0xDB6024D35C80C410uL)),
                ((+1, -12, 0xD237F15597C5D55EuL, 0xA17D1E352E6D1B3DuL), (+1, -5, 0xA0E36E542D7F778CuL, 0x3B84505FC0DD5B51uL)),
                ((+1, -16, 0xD3993B367331C1E7uL, 0x20EC0788CDC428E7uL), (+1, -8, 0xC2D3D851C27969DCuL, 0xC1DE4BE35EB901D6uL)),
                ((+1, -21, 0xB8530FB3981CCD85uL, 0x552FF6773FF3EAD2uL), (+1, -11, 0x9E93104F6792D041uL, 0xDB52B5767A115C25uL)),
                ((-1, -32, 0xC1D3DCDEE9553EC8uL, 0xD0EB17D9DB1696A2uL), (+1, -15, 0x9C2E6D2F9D0C62EEuL, 0x50FB9B803F0687BFuL)),
                ((+1, -37, 0x879988A8271D6505uL, 0xB89F005A8C4DC4D7uL), (+1, -20, 0x8C97F6BD9CF8CA91uL, 0xB68C8B3C916FC906uL)),
                ((-1, -44, 0xFC910DCD5618CD44uL, 0xD4A63AD542C757D9uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xC3A23499C07F5D39uL, 0x58A2287A47C5D1A9uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x9DFDCE1EB6349C30uL, 0xD9DF4DE705AD28FEuL), (+1, 1, 0x88957CF05EE1508DuL, 0x53E52B8B618D0674uL)),
                ((+1, -4, 0xDEB456D36BAB9E55uL, 0x3F7935DD9028FB52uL), (+1, 1, 0x8554046ADC88BE87uL, 0x8D2613106381D22CuL)),
                ((+1, -5, 0xB498B73508CB9BDCuL, 0x93216B7CBF9A77A7uL), (+1, 0, 0x9D35C01440605EDEuL, 0x9070D15228DCD99EuL)),
                ((+1, -7, 0xBA0FC786ECF2FAA2uL, 0x43622B8B950DC1D4uL), (+1, -2, 0xF88F491B0CC57CB6uL, 0x7DD8366710A02686uL)),
                ((+1, -10, 0xFC2436BC85567D0EuL, 0x02E5186919A06A4EuL), (+1, -3, 0x8A1180EB1E5DDBA7uL, 0x280684C428271287uL)),
                ((+1, -13, 0xDF4A6BAE94CF3874uL, 0xDABCA71ABE8C8441uL), (+1, -6, 0xDB4C2107F67D3C73uL, 0x968C6555DC669466uL)),
                ((+1, -17, 0xF607D4085A620362uL, 0x02B620010C06F4C7uL), (+1, -9, 0xF7D75578006887D3uL, 0xBEDC39156B5C5545uL)),
                ((+1, -21, 0x9553E8114D2E7369uL, 0xA52822687861BCB1uL), (+1, -12, 0xC1E2CB472A0E6EF1uL, 0x96D581121D443B65uL)),
                ((+1, -27, 0x914C938DF62A6A17uL, 0x8FE6993458F2CA13uL), (+1, -16, 0xC58E09EBE61C4ED3uL, 0x144B211590C0ED4BuL)),
                ((-1, -40, 0x969163B76A33DEC8uL, 0x06F5B122211F4112uL), (+1, -21, 0xE6DE3BFB892547FBuL, 0x1AABA1A2B1F2B552uL)),
                ((+1, -47, 0xC97A5B7C2E6DE821uL, 0x4B998853654708C1uL), (+1, -27, 0xE0D5F8AC7A15CFE4uL, 0x6D7901C16DE7D871uL)),
                ((-1, -54, 0xB897CA793927585BuL, 0x7674500274AEDBA5uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0x9D23F86327649ADEuL, 0x04EDB643104DEF14uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xA12C882D3078150CuL, 0x89F2491343CE0151uL), (+1, 0, 0xB5C30B273D04A35EuL, 0x7ECE9183171CA079uL)),
                ((+1, -6, 0x928F9FBD64D95151uL, 0x0E7C70E32D076B2CuL), (+1, -1, 0xEB9527CECF75F0EBuL, 0xA5BE3E143A8D5008uL)),
                ((+1, -8, 0x9ADB431094524A56uL, 0x196FC0F8A6483348uL), (+1, -2, 0xB7962E9B854D0515uL, 0x6330C7E51347FF88uL)),
                ((+1, -11, 0xD1153725063CBA57uL, 0xFEC99A22CFEB612BuL), (+1, -4, 0xBEC7840A96236859uL, 0x67AC9FAA85ABE05CuL)),
                ((+1, -14, 0xBAA257AD61913158uL, 0xE8945C19FC7932F2uL), (+1, -6, 0x8A82F5FBB5B2F96CuL, 0x28D593AADDCE41EDuL)),
                ((+1, -18, 0xDC144169ECC0028FuL, 0x671FD035A95973DDuL), (+1, -9, 0x8F186AFBDE514A4CuL, 0x2C854F1A30A03690uL)),
                ((+1, -22, 0xA614AEFE0086C0BFuL, 0xAABEFAB764463F8BuL), (+1, -13, 0xD2031C712018D880uL, 0x802848E66CBD558CuL)),
                ((+1, -27, 0x95825A2BF33D68C2uL, 0x9CA211AC688568B2uL), (+1, -17, 0xD692E087A91F9462uL, 0x2EAA23338271359AuL)),
                ((+1, -33, 0x8BC64017344A16B3uL, 0x6E712924101FA32AuL), (+1, -21, 0x923D597114E87E61uL, 0x960B0596614418BBuL)),
                ((+1, -41, 0xC3CEF9F3E3F19FC5uL, 0x78AE68212FBEB484uL), (+1, -27, 0xF621EE6E736FE97DuL, 0x82DA40224EAB9036uL)),
                ((-1, -57, 0x9A1CA7D2CC294E19uL, 0x7DFD14C3A0EEF989uL), (+1, -33, 0xDD7AF87FA572BCF5uL, 0x7939C38B7F8189EEuL)),
                ((+1, -66, 0xB30A8A5EC940D445uL, 0xC76EADA93CF78600uL), (+1, -40, 0x99076E470CD1D590uL, 0xC458E60583C456C1uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -7, 0xB8937CAB38CD232AuL, 0x1FEC1E80EEA1B3C2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -8, 0xDC52BF6BFDE8FAC9uL, 0xF97DE94AA1F63A49uL), (+1, -1, 0xD61357ED232ACA32uL, 0x873FF70B45379F33uL)),
                ((+1, -10, 0xE5AB4C7579EE8E31uL, 0x8F48E088C18106B0uL), (+1, -2, 0xA105D5F2B60E4EACuL, 0xEEBE88C377763AB7uL)),
                ((+1, -12, 0x88BADA3AB0CBD3EBuL, 0x725FBC7ADBDD5654uL), (+1, -4, 0x8F77413614361407uL, 0x8EAC46058C31A41AuL)),
                ((+1, -16, 0xCC263CD062055081uL, 0x9F98474708C922DDuL), (+1, -7, 0xA7C929FA2572DE06uL, 0x8C086CADC92918EAuL)),
                ((+1, -20, 0xC57FFA664788C538uL, 0x025B2A12D76B1022uL), (+1, -10, 0x86D55ADD1B126890uL, 0x1D76AF8ABB544E4AuL)),
                ((+1, -25, 0xF7390D07166A80C2uL, 0xE24DF67A78F4504BuL), (+1, -14, 0x977BE2386AECE83FuL, 0xFFF69176EA0535F8uL)),
                ((+1, -30, 0xC21B46CEC847E728uL, 0x903A49CEA7F25971uL), (+1, -19, 0xED53032F880CBA8FuL, 0xA6D51A5A92EBFAD2uL)),
                ((+1, -36, 0xB2A486CDAD82576EuL, 0xE71843274EC4715FuL), (+1, -24, 0xFDE96B9A30C95B8AuL, 0x89A2B81A1AD1A5F2uL)),
                ((+1, -43, 0xA8866DBF24F21505uL, 0x62BDA734F8C74194uL), (+1, -29, 0xB1CBBC514735923CuL, 0x072E768122C8326DuL)),
                ((+1, -52, 0xEC7EAE043A9222D5uL, 0x7162A793A39313DBuL), (+1, -35, 0x971AE7E5B6DC4313uL, 0xF04EA644E17A9504uL)),
                ((-1, -70, 0xB91E6E89D2A35687uL, 0xE343FC1E1D201BCFuL), (+1, -42, 0x8784DEECEC9A1DE2uL, 0x08B8C6B70558C6CAuL)),
                ((+1, -80, 0xD6873381CDA96D8AuL, 0xA6371746BE289493uL), (+1, -51, 0xB948CF3803AC5877uL, 0x9E0645FFB021B9EEuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -9, 0xBA05B4CFEDC8049DuL, 0x74AA996F1CF38CA9uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -11, 0xE530673C419355F3uL, 0x756DC90FB90A18B9uL), (+1, -2, 0xDED873EDDCAB7A68uL, 0xDBFD4A1D2A3DCEFBuL)),
                ((+1, -14, 0xF247EC5DD8730D01uL, 0x894EC2142274DC1DuL), (+1, -4, 0xAC58D3492338030DuL, 0xBE8D90027CEE516BuL)),
                ((+1, -17, 0x8F767CFD62C2C437uL, 0x6742261849535222uL), (+1, -7, 0x9BCA07B3FCFAD3FBuL, 0x147C963FF634CBFDuL)),
                ((+1, -22, 0xD08BA81633377861uL, 0x37B4E89AD9435509uL), (+1, -11, 0xB627FAFAC0B7595DuL, 0x6FCE0D3F3A211D1EuL)),
                ((+1, -27, 0xBFD5094723BC15AAuL, 0xD964BDB932CCFB93uL), (+1, -15, 0x8FFB8ACFAC9578BCuL, 0x900109AC1AAE1B61uL)),
                ((+1, -33, 0xDE5CA1CF246C9F3FuL, 0xB6D79D32C415B0A4uL), (+1, -20, 0x9C39C629E6B5AE8AuL, 0x3F3EE78A1A01B326uL)),
                ((+1, -39, 0x9CDAD299A6E34630uL, 0x05339CEC0D443044uL), (+1, -26, 0xE780A7550DFD9E90uL, 0x3F80DB9379981D0FuL)),
                ((+1, -47, 0xFA2AD50BB3BB2703uL, 0xB44356D7337766DAuL), (+1, -32, 0xE4B3A91449844520uL, 0x270A3A0F1F38C835uL)),
                ((+1, -55, 0xC31BA48AE64560B4uL, 0x25CDFD0364FB440FuL), (+1, -38, 0x8FB7A7EB15292B08uL, 0xA0BB641354ED22A6uL)),
                ((+1, -65, 0xD361C445201CE11DuL, 0xC6CD76B915AAA5CEuL), (+1, -46, 0xD3AABD937DAE1B3CuL, 0x7DFDE6ABCC435675uL)),
                ((-1, -87, 0xE5BB77AD33F25AA0uL, 0x936F8459CD3DBC78uL), (+1, -54, 0x9D099350A0F08E68uL, 0x369482FBAAFF540BuL)),
                (Zero, (+1, -64, 0xA5E1157DC1F2B24BuL, 0x816260E7CEDDDF6AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -11, 0xB3C51A79CC925298uL, 0x6B000B60F4CCCDB2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -14, 0xD0C7870F7C934F53uL, 0xFAB0B16D308F6A41uL), (+1, -3, 0xD6617D34BB1D0FD9uL, 0xA22C5A41EB13EF00uL)),
                ((+1, -18, 0xCE178A1C09B5B6C8uL, 0x0B21B8C3D2AB2070uL), (+1, -6, 0x9E574BC5B4893FD4uL, 0xDD11EF18AC5A1F3CuL)),
                ((+1, -23, 0xE16453BFE8EE4E07uL, 0x9BEE6A5CC64920EAuL), (+1, -10, 0x8791395E695CFC57uL, 0xAF229EAC3D0E17DBuL)),
                ((+1, -28, 0x953C1D06C98C75E9uL, 0xE30F9EB8010C0472uL), (+1, -15, 0x94B7752FDE99C72EuL, 0xC856560627F3A0D4uL)),
                ((+1, -35, 0xF5A0B34C8D9A619AuL, 0x05AE30591C4902F7uL), (+1, -21, 0xDA16B8BAAE9E45ABuL, 0x225793655370E812uL)),
                ((+1, -42, 0xF80FECDBE69CF24CuL, 0x7699CBAFC7B56397uL), (+1, -27, 0xD86DA78DC299098DuL, 0x131BA1018AC95E45uL)),
                ((+1, -49, 0x91D3D2589884654DuL, 0x18DDE139083A95ABuL), (+1, -33, 0x8FEDB209BF00447FuL, 0xB32A88A3C323DF92uL)),
                ((+1, -58, 0xB144FA17D8EC7DB5uL, 0x52B8C3A5E45D1D47uL), (+1, -41, 0xF8521F8B84CBA83BuL, 0x87A091957499E54EuL)),
                ((+1, -68, 0xA4CDE2B9A3517FDBuL, 0x4EFD6DD85A4519DFuL), (+1, -48, 0x8234933C302B6B43uL, 0x551878559F24637FuL)),
                ((-1, -89, 0x8EF1AB7120D33DF0uL, 0x24CCC92C72348EDCuL), (+1, -57, 0x923825211F6F9DC5uL, 0x35F5E115FA35D863uL)),
                ((+1, -100, 0xACFA7ECBC5AC1222uL, 0xCF4E5585F6749C5FuL), (+1, -67, 0x8153069ADF50983FuL, 0x3CADC2E9CBE0A54CuL)),
                ((-1, -111, 0x9315C2C1D63B5B36uL, 0x2D8A208BBB2FB5FEuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xADB99CF046DB313CuL, 0xD14F686EE3905524uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xB62EC726F67AEF97uL, 0xCAEBA8CE084C0504uL), (+1, 0, 0x8A0C5B51C42B83D2uL, 0x0649617F42700860uL)),
                ((+1, -2, 0xB3F17ACD8C9C0F66uL, 0x33A7330725B7E505uL), (+1, -1, 0x8B69EB1E0A0E101CuL, 0xB9D56119251A8FA9uL)),
                ((+1, -4, 0xD8A7453D113F52FCuL, 0x85626C6C0C170F19uL), (+1, -3, 0xAA2A087E39107501uL, 0xFD76590C69079EADuL)),
                ((+1, -6, 0xAFA3536734B89D28uL, 0x223DFD944AFB0659uL), (+1, -5, 0x8A8C4BE227657130uL, 0x5DE598623390C025uL)),
                ((+1, -9, 0xC96E551F4406EA68uL, 0x8E66E7A7D9821CD1uL), (+1, -8, 0x9E493F346312424CuL, 0x3DD2E6C6AD32AFBBuL)),
                ((+1, -12, 0xA84A63DB4F844BF9uL, 0x16F9964396E9518FuL), (+1, -11, 0x83B3256D9E97743CuL, 0xD25125EAFFD16231uL)),
                ((+1, -16, 0xD0558A3DD80D79D6uL, 0xC06132D8D7EA4732uL), (+1, -15, 0xA3B288E142F5D41EuL, 0xB7818B4ADD756690uL)),
                ((+1, -20, 0xBFA41744131EB5D6uL, 0x1A907880D3BC01E4uL), (+1, -19, 0x97C9CB59FC6E4F9EuL, 0xC59ED3E10311CA6AuL)),
                ((+1, -24, 0x806E835FBF8CCC9BuL, 0x31039262002AB75FuL), (+1, -24, 0xC5811A52247D87B1uL, 0x2C96DC4523EA2E7BuL)),
                ((+1, -30, 0xE51F0BB22A03C23CuL, 0x3FD801B8F8EB27EDuL), (+1, -29, 0xBA90383B5D90676BuL, 0xFA782D8ACC5D7F0DuL)),
                ((+1, -36, 0xCCA31F42C684476EuL, 0x28FB1ED4CEC4BA5DuL), (+1, -35, 0x95CDD8B9B44945DCuL, 0xB6A6896D00F15A46uL)),
                ((-1, -46, 0x9F49138A161018D6uL, 0xC597169362A04FA4uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA6C4AC2F45683E82uL, 0xA5BC4A9C7E86BDB6uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x989DDF0B96C73338uL, 0x28AE1CBB97501474uL), (+1, -1, 0xED78F873B564007BuL, 0x11232F85FF3A9CC4uL)),
                ((+1, -2, 0x832F77CE548BF441uL, 0x3552E2A02AFE5832uL), (+1, -2, 0xCDACA72E234F599CuL, 0x998E051281BFF9E0uL)),
                ((+1, -4, 0x8CFB2AC11D4DF6F6uL, 0x311037136512714EuL), (+1, -4, 0xDD96940F182005CEuL, 0x4851AD3695907A91uL)),
                ((+1, -7, 0xD4DD468BD1D82196uL, 0x91016F372727C8EDuL), (+1, -6, 0xA7416ECF93FB8E6BuL, 0xF8BE24E268F22B1EuL)),
                ((+1, -10, 0xF0848D8B26DCE889uL, 0x44F0C479273C437FuL), (+1, -9, 0xBCE0374200C25FD1uL, 0x974712DD2F31E01CuL)),
                ((+1, -13, 0xD361899DA8C8A337uL, 0x5E60D5EBD1A08195uL), (+1, -12, 0xA5FF03F04927CF75uL, 0x763416A88360485CuL)),
                ((+1, -16, 0x9426144B6FCF51B2uL, 0x2F285C9143A1EBCEuL), (+1, -16, 0xE8B9C03351A4EF63uL, 0x3118947357280CEAuL)),
                ((+1, -20, 0xA84A5589F303DC24uL, 0x85FF51D27A3EA7C5uL), (+1, -19, 0x842F2BDA8637B641uL, 0x1A0BA8F01F1F0014uL)),
                ((+1, -24, 0x9C5BB087A48D781CuL, 0xA8EE28B3C21859E2uL), (+1, -24, 0xF59AE1CCFC666EE1uL, 0x639EEF9087F17A44uL)),
                ((+1, -29, 0xEE527E06400D77C2uL, 0xF0D7DF047D2A9255uL), (+1, -28, 0xBB20297FD81D71FBuL, 0xD0EAF6DDF816A8C9uL)),
                ((+1, -33, 0x9458E3D2D3A4EDC7uL, 0x0156AC2667C384E8uL), (+1, -33, 0xE92DB5BF681EEC88uL, 0xEA6DB28D2884C197uL)),
                ((+1, -38, 0x94BDBCDC6BD994A7uL, 0x803362920D1A8E08uL), (+1, -38, 0xE9685949DE4AEE55uL, 0x933CDE51D502A629uL)),
                ((+1, -44, 0xE8B8DE7C57802600uL, 0x4A916409EDBF1E47uL), (+1, -43, 0xB6FF74D682DCE7E2uL, 0x1FDD106CE001C182uL)),
                ((+1, -49, 0x83FAA694F202365EuL, 0x84DCE11E3ABB8B00uL), (+1, -49, 0xCF0E2FC267403B47uL, 0x733AE6DE062F6F94uL)),
                ((+1, -56, 0xBB0D01FCDF3D0834uL, 0x7222E9095B5923A0uL), (+1, -55, 0x931654FAB1FBBB2BuL, 0x84B4E45587722FCCuL)),
                ((+1, -73, 0x91FBC5BF2E871FD3uL, 0x44785B6AF06287B6uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA301C18B3FBABE4BuL, 0x6736743D5D63299EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x85C1376DAFD16CC5uL, 0x6D318C19883EBC78uL), (+1, -2, 0xD21F8927ED383E32uL, 0x7E860F04B81DD7AEuL)),
                ((+1, -5, 0xD590C17FFADB0C9EuL, 0xED0C0C3998B08FDDuL), (+1, -4, 0xA7B9F936F2C4E7EEuL, 0x12AE41C8D108B4A3uL)),
                ((+1, -8, 0xDCCE92881E975136uL, 0xA2C8D56ED83601AAuL), (+1, -7, 0xAD6C50DB52A723F9uL, 0xD6D005221406BC02uL)),
                ((+1, -11, 0xA5DD91AD171BA9B8uL, 0xCA18779C705293FDuL), (+1, -10, 0x82457E64D636F3D4uL, 0xF020C58E084E575DuL)),
                ((+1, -15, 0xC06A11F8C07F9E07uL, 0x9DEE182B48D0AB11uL), (+1, -14, 0x971E9EDFE5AE6ABEuL, 0xFBFABC957105E3C9uL)),
                ((+1, -19, 0xB2B71CBF657BFF67uL, 0xDF4A0F8736FA75FAuL), (+1, -18, 0x8C5D6A418E8F980DuL, 0x1E4901B4F3DA01CEuL)),
                ((+1, -23, 0x87E7F9CDED0FAA5DuL, 0xCDD3D49761323D52uL), (+1, -23, 0xD57A584C26EFE37AuL, 0xF253515E580A03AEuL)),
                ((+1, -28, 0xAB79930168773708uL, 0x57B5F42D7326D2FCuL), (+1, -27, 0x86AD718AAEA8C42EuL, 0x3E5B4DEA030513B1uL)),
                ((+1, -33, 0xB48B074CEA40993DuL, 0x1B2D54D2BB16D096uL), (+1, -32, 0x8DCC11FDB6FFA926uL, 0x1E0CFE0BA09B556FuL)),
                ((+1, -38, 0x9E7C0BB3EAE14036uL, 0x5AB25F2BD1672E66uL), (+1, -38, 0xF8F2B844C8E0BDB0uL, 0x192FC4888F33679BuL)),
                ((+1, -44, 0xE5E5A60CDF341E9AuL, 0x7EBA78282455B3A5uL), (+1, -43, 0xB48F6C407E619BE0uL, 0xCF707170D7DDEB1EuL)),
                ((+1, -49, 0x8708DEE7EC372795uL, 0x203B7315F0E6769EuL), (+1, -49, 0xD41CBB0E78B9E1FFuL, 0xD00C2C370C4723EAuL)),
                ((+1, -56, 0xF695BC6D30CCB392uL, 0x8B5B0073F3E5A8F0uL), (+1, -55, 0xC1AACE4C5B6FCC1AuL, 0x1858F66EEF597902uL)),
                ((+1, -62, 0xA0FFC34BFD58E29EuL, 0x74942036C7C21633uL), (+1, -62, 0xFCE595AA72B2AD73uL, 0x1D0D9C97008B045DuL)),
                ((+1, -70, 0xEDF3A123B7124D71uL, 0xFB35DF6F71CD9372uL), (+1, -69, 0xBAE3051DC25C71BAuL, 0x7DAF6D4F7A87BB76uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9837F89689087uL, 0x763B33B78FEBBCC2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF46D3FFC82BB0FE8uL, 0x59246A328C8C3CB1uL), (+1, -2, 0xBFF8EA734DBC92E2uL, 0x483CB336AB417FC0uL)),
                ((+1, -5, 0xB41CCC04B90A11DCuL, 0xF594AE00E911743CuL), (+1, -4, 0x8D75C3A9CD7A3D61uL, 0x7381E27990532CC0uL)),
                ((+1, -8, 0xAD598A31E96DE2E1uL, 0x9EC1E7BC86EC4928uL), (+1, -7, 0x882609D230DD0E3EuL, 0x0DB70EB296E4702CuL)),
                ((+1, -12, 0xF459D0BDA41C055EuL, 0x3897004C5877D0CBuL), (+1, -11, 0xBFE9A6D19119C4B1uL, 0xAEC73714E825FC21uL)),
                ((+1, -15, 0x8607DDBB6DDA88ACuL, 0x23ACEDAF5030E523uL), (+1, -15, 0xD288F40FCAABB771uL, 0xD54025E253E55B67uL)),
                ((+1, -20, 0xED690133262CE024uL, 0xBB7F885F7F38D710uL), (+1, -19, 0xBA7625D77F048AB8uL, 0x5D85D1670D686689uL)),
                ((+1, -24, 0xADBD13160CCDA68EuL, 0x9A39533DD208B554uL), (+1, -23, 0x88743672A1EE4075uL, 0xF494C17DB3281487uL)),
                ((+1, -29, 0xD53472F22CC18996uL, 0x1F3CF27E46963204uL), (+1, -28, 0xA7736267B5604906uL, 0xFBF298EA065F6417uL)),
                ((+1, -34, 0xDD2D24578081C1EEuL, 0x397B8BE73C7CD98BuL), (+1, -33, 0xADB62412DDE4FE82uL, 0xF17CAF1C9285D8C8uL)),
                ((+1, -39, 0xC27BD6FF25DB0405uL, 0x920AE887141E52FEuL), (+1, -38, 0x98BF47284E995FE6uL, 0x759ACB1808A1F179uL)),
                ((+1, -44, 0x9083A2580B600E53uL, 0x6BF414DEFEA94EB1uL), (+1, -44, 0xE3009B32B60733FAuL, 0x0D6005EE3A3304CFuL)),
                ((+1, -50, 0xB3A8182B6CE4FD01uL, 0xBD206DF564971AC4uL), (+1, -49, 0x8D1A1B4260C69B59uL, 0xD185FDA485524E28uL)),
                ((+1, -56, 0xB6FE96E1BE431FC6uL, 0x9E6BE246F364F489uL), (+1, -55, 0x8FB939AACDB2FAE8uL, 0x462E03093ABA5103uL)),
                ((+1, -62, 0x92C5029B88CAE5C6uL, 0x2987D949D776D1C1uL), (+1, -62, 0xE68B8BDFF6E39456uL, 0x15B591E6ABFF7DD9uL)),
                ((+1, -69, 0xAA9DD89F131B2C13uL, 0x42832F941C197FA1uL), (+1, -68, 0x8600800B10CD7725uL, 0x145647156D7C8235uL)),
                ((+1, -77, 0xE79EB502FC39D547uL, 0x807AC6B7605B6AB7uL), (+1, -76, 0xB5E9F43BC9ADF48FuL, 0xFC67B7CB145D2EBDuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E44154DuL, 0x31EBF782DC89ECD4uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8EA2B5E57A8A456FuL, 0x67161DE2C166A688uL), (+1, -2, 0xE00D2C55744537FCuL, 0x85E618A34CC3E2C5uL)),
                ((+1, -5, 0xF5D5861C489EB50FuL, 0x542404D9BC23A8A2uL), (+1, -4, 0xC113DFC29614F0ADuL, 0x5E617E7A771D23BBuL)),
                ((+1, -7, 0x8AD6C88B0F1780C1uL, 0xAEF02B521E72FB15uL), (+1, -7, 0xDA16793B879E5388uL, 0xB76548E4C37462B1uL)),
                ((+1, -11, 0xE6D0AC9D4C7C9FABuL, 0xE3187D93204DF022uL), (+1, -10, 0xB54822E17F6FDD75uL, 0x42D8DF70F7D4C049uL)),
                ((+1, -14, 0x9645C43FEDF81C7BuL, 0x5066E6CBFD88951CuL), (+1, -14, 0xEC0C2B068CEFBB49uL, 0xEEA4815D2CAD2884uL)),
                ((+1, -18, 0x9F38366D25CB4B01uL, 0xA3A9622D042D9D55uL), (+1, -18, 0xFA19FE032F893D63uL, 0xEA8439C7E0E51030uL)),
                ((+1, -22, 0x8CBED1248EAC0544uL, 0x06F6066270C74568uL), (+1, -22, 0xDD15132F2CC2302FuL, 0xD0BE32CB4BBE5909uL)),
                ((+1, -27, 0xD309D21E93E16457uL, 0xD1F3FD2BF54896FAuL), (+1, -26, 0xA5BFC7C956964371uL, 0x8677BA3575694EE5uL)),
                ((+1, -31, 0x87C705C671C92079uL, 0xBD6190F171A544DFuL), (+1, -31, 0xD547584F3D18E2B7uL, 0x07FA5FEFB139E857uL)),
                ((+1, -36, 0x95FA795E7F21F8CEuL, 0x3790489D5B64D891uL), (+1, -36, 0xEB95E61D2B921D28uL, 0x2973C1F0CDFB19DBuL)),
                ((+1, -41, 0x91B82B32FE5CB522uL, 0x7A5170F7604ADD41uL), (+1, -41, 0xE4E540516019A6EAuL, 0xA3E926BBCBCC8340uL)),
                ((+1, -47, 0xE31995739DCA52F7uL, 0x7A37E15AB93180F7uL), (+1, -46, 0xB25D26CB25C25C78uL, 0x2158BB1E381FBE61uL)),
                ((+1, -52, 0xB8AC0BDD20FCF0E9uL, 0x2AC4F0D0E490DD34uL), (+1, -51, 0x910A851CC7EF1784uL, 0x3A7214DF81E0F182uL)),
                ((+1, -58, 0x86D7E7935F43AE36uL, 0x4145BD248991074CuL), (+1, -58, 0xD3CFBD4511B78BFAuL, 0x6B809687E7F37FD8uL)),
                ((+1, -63, 0x87692C8C99443F33uL, 0xA26B576E4A574C20uL), (+1, -63, 0xD4B3ED92A57AE89CuL, 0x5F47CF6A155017CFuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_0 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xE315A1E749CD06E3uL, 0x4A2444091D003670uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0x8C98979686DC559AuL, 0x549CEC83030D5C38uL), (+1, 1, 0xFBB37FE6CB7650B1uL, 0x9183B43BE5792C62uL)),
                ((+1, 1, 0x9D242795B2C8EDC9uL, 0x7B0F84182DF76E4EuL), (+1, 2, 0xFD19135F0A129527uL, 0x299BA28BCFF7E15FuL)),
                ((+1, 1, 0xD1C4DD3487EBA2BBuL, 0xD2805CEB0E39E12FuL), (+1, 3, 0xA8A8A7A6FF6EEC9EuL, 0xE3D3125172197F78uL)),
                ((+1, 1, 0xBAC41BE00D894F50uL, 0xC5ECC921C5C0A3CEuL), (+1, 3, 0xA505522A4DA26C1BuL, 0xB32D6DBEEBA8DDE2uL)),
                ((+1, 0, 0xEB165C57AD416536uL, 0x7EE5BC0989AC615AuL), (+1, 2, 0xF981F172E1A2FB25uL, 0x6265490808C8D8F6uL)),
                ((+1, -1, 0xD8A13BC8599E7DA5uL, 0xF5327410A25FF8DFuL), (+1, 2, 0x95EE538AEE8B67B0uL, 0x7671C0B88125DAFCuL)),
                ((+1, -2, 0x95319684ACBF1E16uL, 0xF744AA7DA0E84BFCuL), (+1, 1, 0x915BBCC266B65814uL, 0x301B11C3EE332AF3uL)),
                ((+1, -4, 0x9A71CA93351F1803uL, 0x8A1E34DB76EE8FEFuL), (+1, -1, 0xE4982BB0E35F0383uL, 0xB1CFAD267BFE6154uL)),
                ((+1, -7, 0xEC9EB78D31625C13uL, 0x267DB36036E08F3BuL), (+1, -2, 0x915399F5FE8BDBB9uL, 0x67ED7EBE28735840uL)),
                ((+1, -9, 0x805295ADE55DEA57uL, 0x192585B8C931B427uL), (+1, -4, 0x938E3A9FFAD6B347uL, 0x16BF23E8A1D90A8FuL)),
                ((+1, -13, 0xB2D4D7CA61100700uL, 0x96951CF68693323DuL), (+1, -7, 0xE93BCF7A48666583uL, 0x15EB5893F76C0F15uL)),
                ((+1, -18, 0xED5CAFFFE7E10F74uL, 0x0318A1A42E3D8CE7uL), (+1, -9, 0x88CB85476147F95AuL, 0x3EA6F7FC74F2A0C2uL)),
                ((-1, -28, 0xD2AE004F6F63DBBBuL, 0xD007A474FAC8FEBCuL), (+1, -13, 0xD809CC0C1278CA07uL, 0x50F0EC56C3F3D064uL)),
                ((+1, -34, 0xE95529BFF39D51AEuL, 0xA7FF6153176E193EuL), (+1, -17, 0xB2CF035EB4490447uL, 0x65D4D88B6CA1837FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1p5_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xACF73FADBECAA9A5uL, 0x20F452296A8D4861uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x8250D6A20C820A12uL, 0x1C4C43D025E72E6CuL), (+1, 1, 0xBA52C4B8F5A183A6uL, 0x3D28AD95D2AC864AuL)),
                ((+1, 0, 0xADC9DDB13A771C6CuL, 0x5E2A96E9C610BAABuL), (+1, 2, 0xA5EEB648FFD73103uL, 0xFA3C60869F8772A7uL)),
                ((+1, 1, 0x8685B7390D90C68BuL, 0x323F7F012F71D7FCuL), (+1, 2, 0xC80AA8AE9D37DA7DuL, 0x76732CFCACC25A7BuL)),
                ((+1, 1, 0x855DA27EBE7B4EDFuL, 0x8BF1D473CD56601BuL), (+1, 2, 0xB9C69A55B2FE30D7uL, 0x7A09DA58673386D0uL)),
                ((+1, 0, 0xB0B2FEC02FED53D7uL, 0x818D863A8ADF915EuL), (+1, 2, 0x865894E466EF09DBuL, 0xE946297AE54C9A45uL)),
                ((+1, -1, 0x9F6F5DDE83CEE657uL, 0x60038BDE6318C2BCuL), (+1, 1, 0x9E3E9F199414DAF7uL, 0x9CBFA5BFE06E2495uL)),
                ((+1, -3, 0xC67915F2F75CAA6DuL, 0x21B4B3A5E0550077uL), (+1, 0, 0x958622EE982F906AuL, 0xB734D67EA6135BCBuL)),
                ((+1, -5, 0xAC8CBA27BC560B7CuL, 0x378F1A2E29386B5EuL), (+1, -2, 0xE81C649F6DDA529DuL, 0x821011A9AC48FA81uL)),
                ((+1, -8, 0xC8BB7F073D760DC0uL, 0x790F70CB23E8A1A9uL), (+1, -3, 0x8D04C775ED85E794uL, 0xD04E740F2743CF0EuL)),
                ((+1, -12, 0xE2D41D2F7D3E767AuL, 0x54A19EA712268944uL), (+1, -5, 0x88D94189D780474AuL, 0x9485E9204DBC01A2uL)),
                ((-1, -25, 0x9A63DEA751AEF1ABuL, 0x7EADAEB23821B362uL), (+1, -8, 0xB4686B9C0E0EA7C9uL, 0x52B1563343018344uL)),
                (Zero, (+1, -11, 0xAAA8B7CD3F9C2D31uL, 0x7A99DD654F761256uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1p75_1p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0xFF6F93D64B8F16FDuL, 0x50CB5CA841EAFDB0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xE0E4D1E4D4E3F221uL, 0x83982C6E41753F56uL), (+1, 1, 0x843A3BE5E03884FCuL, 0x5DB7D32EE500AEB9uL)),
                ((+1, -1, 0xAD7254E93536FACBuL, 0xC943B54336DE4BBFuL), (+1, 1, 0xE40791D1E567FFF7uL, 0x611B57D2D1E440B6uL)),
                ((+1, 0, 0x98671CB3F91F6EDEuL, 0x26BA7AADC694A2B2uL), (+1, 1, 0xEADB1D71D1FE2F23uL, 0x5B28A7EEEF620FC4uL)),
                ((+1, 0, 0xA59F85632C9B3327uL, 0xE742ABD983B2AB90uL), (+1, 1, 0xD67989090080809BuL, 0xA9E2B10E6E3057B0uL)),
                ((+1, -1, 0xE0AE3CE18D470F4AuL, 0x9342BCED97C4CE78uL), (+1, 1, 0x88097D44A146BA83uL, 0x2E0FDDA68E29895CuL)),
                ((+1, -2, 0xB3EAAA8885F6016BuL, 0x28360DE1BE5E920EuL), (+1, 0, 0x9FBA6AB860D72043uL, 0x2654C5A369F9857FuL)),
                ((+1, -4, 0x8E0C7555D63CCBF6uL, 0xC508C52AEE26A6AAuL), (+1, -2, 0xFD3940FF27B8623BuL, 0xBE2A1FD5E37ED142uL)),
                ((+1, -9, 0xDC0B3FA929FA5A67uL, 0xC95815333F3C8868uL), (+1, -3, 0xC7241EECFED0B051uL, 0x2C062705BCF78025uL)),
                ((+1, -12, 0xB57BAF2D5FCFC0E8uL, 0xD78AF270367C436BuL), (+1, -5, 0xA43886B060CA276BuL, 0x071B9B7FC657B13FuL)),
                ((+1, -13, 0xCC75B75868C3ABE4uL, 0xE854D7E5635C3A3BuL), (+1, -7, 0xA9572EA02DA643A6uL, 0x410273CD4830E921uL)),
                ((-1, -15, 0xABC4422B4E918DF4uL, 0x33AD46E3773F1671uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_1p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -8, 0xD53E03A99B7FE181uL, 0x7B6797698939ED63uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x81874523EF0195CAuL, 0x60EDFC7CC9C73B7AuL), (+1, 0, 0xFD8AA46DCF1F40E3uL, 0xD1B42957C3EDC84DuL)),
                ((+1, -2, 0x8D279F5128E7ED0CuL, 0x1801D969CE0AD1E5uL), (+1, 1, 0xF1C8EC181FFC3035uL, 0x4BE29A62AEF0D948uL)),
                ((+1, -1, 0xB4B5AFE5FDC154F8uL, 0x61168AFF92C5766FuL), (+1, 2, 0x83FAAC8618CE45ABuL, 0x3ADF12788D509F6AuL)),
                ((+1, 0, 0x951DD9D182DF6E7FuL, 0xE0C3AB53734969BEuL), (+1, 2, 0x8780D61596CD8D13uL, 0x623135867128E00BuL)),
                ((+1, 0, 0xA30E921291247933uL, 0x810E9C80808A5C1DuL), (+1, 1, 0xC225719B11D586F3uL, 0x300B80F19D600C7EuL)),
                ((+1, -1, 0xE7A3DE6DB196B5E0uL, 0xD2EBEE4664059C01uL), (+1, 1, 0x82559937D1C3BDD5uL, 0x961CEAC22DF248B7uL)),
                ((+1, -2, 0xC34A2168B4564B1AuL, 0x69FDE3470103F8A5uL), (+1, 0, 0x804AF8C98F0FFAA4uL, 0xFA191F7D68C7FC34uL)),
                ((+1, -4, 0x8DAA7DD4C70A0EA8uL, 0x97187DBBDD15C1E5uL), (+1, -2, 0xE53FFDB467BE3012uL, 0x010D31B78B08534EuL)),
                ((-1, -8, 0x97DEECE137C13B2EuL, 0xF537E0265360EB5AuL), (+1, -3, 0x995969BB92EA5800uL, 0x78248B50AAF0DBFBuL)),
                ((-1, -13, 0x83F9E790ED6E5F8BuL, 0xB352C5EB2139A31DuL), (+1, -5, 0x99E47B4452857F77uL, 0xDD09DFA09A4EA1E5uL)),
                ((+1, -10, 0x8C068E3D91D3C5F7uL, 0x85E89B6596D1749DuL), (+1, -7, 0x86AFF7B6AF2060F2uL, 0xA4484A83DA301D85uL)),
                ((-1, -13, 0xC89CC543126E4A73uL, 0x7D3AEB349876452DuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA1917E8D995A6725uL, 0x8D4D53E4F7CDA88BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xCF0774BE265D47E1uL, 0x353E8D2F506DA6E7uL), (+1, -1, 0xA6AD84A1A4EAE260uL, 0x11418D3218BADC68uL)),
                ((+1, -2, 0xDBFF30DA37864859uL, 0xEFDA3ED1C7EA61F8uL), (+1, -1, 0xAE2DB82E3D277E32uL, 0x662DD2AC4B4D98A6uL)),
                ((+1, -3, 0xDE7A60D35F4E65D5uL, 0x8D456CD14DD6B1F6uL), (+1, -2, 0xB2FF80B61E58E75BuL, 0x021EE30078E0346AuL)),
                ((+1, -4, 0xF95D14FA17320EABuL, 0xF99C87F4718642FBuL), (+1, -3, 0xC580B57C9CFE1166uL, 0x3D23B017FF4BEB80uL)),
                ((+1, -5, 0xC476F3CC039773FDuL, 0x28DB58F1DE1DEB69uL), (+1, -4, 0x9DE5D7CD8EE5D9E8uL, 0x8E20A267A13E11A8uL)),
                ((+1, -6, 0x94DD340064CE5601uL, 0x7BCEF32371053D5FuL), (+1, -6, 0xEC0AB37708291327uL, 0x3A66E94F0B61D255uL)),
                ((+1, -8, 0xB45029198C18398AuL, 0x7A904E1D389F6DA1uL), (+1, -7, 0x90B77B76B5B22E9BuL, 0x55C0D7E464A8BE0AuL)),
                ((+1, -10, 0xC52A82E3757AB73FuL, 0x66640F86C6D7D105uL), (+1, -9, 0x9C8592C20782B408uL, 0x8C2414AC8D19F74AuL)),
                ((+1, -12, 0xB2CDF88600D794F4uL, 0xCD0E54ADD0966CDDuL), (+1, -11, 0x8F495B68B0E64729uL, 0x1515419B49A3ADA1uL)),
                ((+1, -14, 0x8C198A666AD10109uL, 0x08BF0F0D2875BC1DuL), (+1, -14, 0xDEC4D322CBEBA4C9uL, 0xB289FF52B1FF480EuL)),
                ((+1, -17, 0xB3C2B88C5D1FA24AuL, 0x51A60A306A65DF5DuL), (+1, -16, 0x8FD4E1D881DA0926uL, 0xF44E05ABCD28AC62uL)),
                ((+1, -20, 0xBBB4F6A4888BD6E8uL, 0xBD112317C512496FuL), (+1, -19, 0x95770B58493661A5uL, 0xDE4BD138808CEFA7uL)),
                ((+1, -23, 0x94F596FD28351367uL, 0xDE7A9FAE70E25215uL), (+1, -23, 0xEE0584992C62D5D3uL, 0x8A53FEB69EA16A7AuL)),
                ((+1, -27, 0xA4D60C1A308C927CuL, 0xC0F64A74006D6B8FuL), (+1, -26, 0x836AF3C3EB778623uL, 0x7F7A2B67B8A4D4A4uL)),
                ((+1, -32, 0xC258D98E84B223EDuL, 0x6BD6A2C21E2B3CE8uL), (+1, -31, 0x9B20E404F03B9C30uL, 0xD42031D681110BBEuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA07A30CEDEEF2742uL, 0x731FD88D98E21CEFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xB74141BD22C5F69BuL, 0x1DE6FBABEBDB4393uL), (+1, 1, 0x9233762C35CC591AuL, 0x9F210F50C08C5194uL)),
                ((+1, 0, 0xC2D0963A329705D9uL, 0x6C3EAF274DC221F6uL), (+1, 1, 0x9B702ABBB73DB64FuL, 0x8CF1E1B02D936D7BuL)),
                ((+1, 0, 0x814E0CA1C66CA9F9uL, 0x33E756B7DAD27B73uL), (+1, 0, 0xCE5791148B43AA76uL, 0x8E1C1C676EE763A0uL)),
                ((+1, -2, 0xEEE0810F477BE149uL, 0x329A0DBFB7ED8F1DuL), (+1, -1, 0xBE98EBF83ACF748FuL, 0x6331DF48EAAC21EBuL)),
                ((+1, -3, 0xA264733F89BF362CuL, 0x8B46BA9686B8A4BDuL), (+1, -2, 0x8191E1492FB02629uL, 0x87E38DFA575BDDBFuL)),
                ((+1, -5, 0xA7DD4CC657B68209uL, 0x41B8BBB2BD4B03DCuL), (+1, -4, 0x85EFB96492E33EFBuL, 0xB3E00C90E2FC0228uL)),
                ((+1, -7, 0x86CB73E3A5DF88FDuL, 0x8F05694021AB4B4EuL), (+1, -7, 0xD719ED75DBD4351CuL, 0x4A864E693E997A28uL)),
                ((+1, -10, 0xAA8AD0915D1668E5uL, 0xF7C5C826AD4A704BuL), (+1, -9, 0x8812E1CBC900EDF2uL, 0xC2E0603302D10BAAuL)),
                ((+1, -13, 0xAB40AE12BE93EFCCuL, 0x7EDF832EC03BA07CuL), (+1, -12, 0x88A3519BBC2583B1uL, 0x181D1112244585CCuL)),
                ((+1, -16, 0x884F76217AB70C33uL, 0xBC0C9D66DCD02EAEuL), (+1, -16, 0xD986616F731947F3uL, 0xE5FEDAED0F16879EuL)),
                ((+1, -20, 0xAA0BCF76F7F5A9DCuL, 0x6DDE0E191A6F357EuL), (+1, -19, 0x87AC431120EB8FE9uL, 0x7405B340ACA0E33EuL)),
                ((+1, -24, 0xA13AD71F5747BEDEuL, 0x79E32D9335E61C36uL), (+1, -23, 0x80A5DB82DBAC3170uL, 0x0AE896B1D38C8D30uL)),
                ((+1, -29, 0xD6E0F2A27906E80CuL, 0x16B825E80E9C8252uL), (+1, -28, 0xAB70A450F5B07514uL, 0xE7FA678203531758uL)),
                ((+1, -34, 0xAB40F49E73692B41uL, 0x3508636BABFC80EAuL), (+1, -33, 0x88A66B8A688D1732uL, 0xEDF5FC7A37E827E2uL)),
                ((+1, -53, 0x83CE2E4BB6CA9808uL, 0xEAB0488448A3BAD2uL), Zero),
                ((-1, -60, 0xA8E4FB17EFE47017uL, 0xCBE0AC1C40784C34uL), Zero),
            }));

            private static ddouble PlusValue(ddouble x) {
                Debug.Assert(x >= 0);

                if (x <= 64d) {
                    ddouble y;

                    if (x <= 1d) {
                        Debug.WriteLine("pade minimum segment passed");

                        y = ApproxUtil.Pade(x, pade_plus_0_1);
                    }
                    else if (x <= 2d) {
                        y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                    }
                    else if (x <= 4d) {
                        y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                    }
                    else if (x <= 8d) {
                        y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                    }
                    else if (x <= 16d) {
                        y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                    }
                    else if (x <= 32d) {
                        y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                    }

                    return y;
                }
                else {
                    int exponent = double.ILogB((double)x);

                    ddouble v;
                    if (exponent < 8) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -6)), pade_plus_expp6_8);
                    }
                    else if (exponent < 16) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -8)), pade_plus_expp8_16);
                    }
                    else if (exponent < 32) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -16)), pade_plus_expp16_32);
                    }
                    else if (exponent < 64) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -32)), pade_plus_expp32_64);
                    }
                    else if (exponent < 128) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -64)), pade_plus_expp64_128);
                    }
                    else {
                        v = 2 * RcpPi;
                    }

                    ddouble y = v * Square(1d / x);

                    return y;
                }
            }

            private static ddouble MinusValue(ddouble x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 2d) {
                    ddouble y;
                    if (x <= 1d) {
                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0);
                    }
                    else if (x <= 1.5d) {
                        y = ApproxUtil.Pade(1.5d - x, pade_minus_1p5_1);
                    }
                    else if (x <= 1.75d) {
                        y = ApproxUtil.Pade(1.75d - x, pade_minus_1p75_1p5);
                    }
                    else {
                        y = ApproxUtil.Pade(2d - x, pade_minus_2_1p75);
                    }

                    return y;
                }
                else if (x <= 8d) {
                    ddouble v;
                    if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }

                    ddouble sigma = Exp(x * pi_half - lambda_bias);

                    ddouble y = v * Sqrt(sigma) * Exp(-sigma);

                    return y;
                }
                else {
                    return 0d;
                }
            }

            public static ddouble Value(ddouble x) {
                return (x >= 0d) ? PlusValue(x) : (x <= 0d) ? MinusValue(x) : NaN;
            }
        }

        private static class CDFPade {
            private static readonly ddouble pi_half = Ldexp(Pi, -1);
            private static readonly ddouble lambda_bias = (+1, 0, 0xB9CD764B54262B88uL, 0xCAD8DC6940262974uL);

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA27FB769C55610C7uL, 0x8A68CD7CD3DB7F2CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xBA660E88ECD7E72AuL, 0x69531847E7CCF51BuL), (+1, 1, 0xAD4427B5D8A42F38uL, 0x0A0F0E8D75B5E628uL)),
                ((+1, 0, 0xC7E8C7E4885579E5uL, 0x886EAD90933F4B1EuL), (+1, 1, 0xE0B91C56D10E51D6uL, 0x8C465650200505B5uL)),
                ((+1, 0, 0x83245E9061376340uL, 0x3524204393DE856AuL), (+1, 1, 0xB6A732F5B94D34E4uL, 0x96656E3F94AC7C37uL)),
                ((+1, -2, 0xE8A69E0330ED37D0uL, 0xFDD82CBD1C91280BuL), (+1, 0, 0xCDA37FE617BFA1E1uL, 0xE70BB55D8DE36562uL)),
                ((+1, -3, 0x917F0C28FB8A94BEuL, 0x59BFCBB87DE9D9D5uL), (+1, -1, 0xA7A3F04821F81A50uL, 0xE6E8547A136439C2uL)),
                ((+1, -5, 0x816C816C3651FD9DuL, 0x56438B35AE727567uL), (+1, -3, 0xC91F6D4BDD056D12uL, 0x43ACEF9A89CB40E2uL)),
                ((+1, -8, 0xA06E5D0C47B73DA5uL, 0xFBE3F414D4C88A48uL), (+1, -5, 0xB0A8336B724C8A3CuL, 0x398AD3FB40A7EA97uL)),
                ((+1, -11, 0x82A1C2F567785274uL, 0xE8A7026AE52ACEE1uL), (+1, -8, 0xDCF4DBFEE8040A02uL, 0x058251E97721DD12uL)),
                ((+1, -16, 0xF388994572F98E25uL, 0x7EB15792EC0C7DB7uL), (+1, -11, 0xB8BBE5123EDC3140uL, 0x1624146D536D11C2uL)),
                ((+1, -21, 0xB30422861A263A01uL, 0x7CB337139ABD762EuL), (+1, -15, 0xB3D7CFEC05F59179uL, 0xC1E08B38383AA808uL)),
                ((-1, -35, 0xFF5DF0516A90099FuL, 0x404702961F337401uL), (+1, -20, 0x8B37E3AB922656F5uL, 0x04AA08093B3628D8uL)),
                ((+1, -42, 0xDC7D294165689959uL, 0x7EDBE07C62262E0AuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD821D91BC8DE6DC7uL, 0xE4E97F1B9F56AD75uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x8B1461CA3B270C29uL, 0x5367EF7B50D7DCDDuL), (+1, 0, 0xD6520C546D11F3A6uL, 0x61F0A293A7DBFA08uL)),
                ((+1, -2, 0x9E5CDFC0E26C1C2FuL, 0x0676680A2DCC7BB5uL), (+1, 0, 0xA363AEBD7BC24772uL, 0x91B96863771289DDuL)),
                ((+1, -4, 0xCA8471DDEE238246uL, 0x6C82E59F14A72C98uL), (+1, -1, 0x91D6BE364754099DuL, 0x48603EF3D29C6D98uL)),
                ((+1, -6, 0x97BB691AC40270F3uL, 0x9F82F0B3940E5DD3uL), (+1, -3, 0xA3C24D1653DBF2ABuL, 0xE3C508008C1F4B81uL)),
                ((+1, -10, 0xEBB2E0E6D63D6AEDuL, 0xF95F9E814CA266D1uL), (+1, -6, 0xE43CA9DFC624DD39uL, 0xE1498EB0839D8900uL)),
                ((+1, -16, 0xCC41AF3F23D796BAuL, 0x832ED9FD8DF0798FuL), (+1, -9, 0xAB12556C0C671A3DuL, 0xBD098EE08F925273uL)),
                ((-1, -17, 0xC41837BF4E0079A7uL, 0x9F413B4093A8A2BBuL), (+1, -15, 0x8E33AD1E49984626uL, 0x98DBAE869F4A174CuL)),
                ((-1, -21, 0xD412CBEB852068B1uL, 0xEF83C8BE23F82D54uL), (-1, -16, 0x8F75350CB59DE7D9uL, 0x17D2B91496C2B2BEuL)),
                ((+1, -32, 0xF703E751B0D1263FuL, 0x47B8F64B8CCC8287uL), (-1, -20, 0xA28EDAD39B8B4600uL, 0x14BCEF52A701F5CDuL)),
                ((-1, -37, 0xCBFA0B50964755F7uL, 0xA51D2FAA698D4319uL), Zero),
                ((+1, -43, 0xD41757AFB6ACA9C4uL, 0xF5E3DD27614F1675uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x977F2C9F53C916FAuL, 0xC7E702D675543F01uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xC45E049EEA227FBDuL, 0x96C5DCB679C2377AuL), (+1, 0, 0xCF3BE856B3F6A376uL, 0xB3007D59DC3322D1uL)),
                ((+1, -3, 0xE3206BD38B1F2EC5uL, 0x665CECA9CB4AA3EDuL), (+1, 0, 0x98214B9B4B84301CuL, 0x795379BD907D59ABuL)),
                ((+1, -4, 0x983937D2210DFAB9uL, 0xFC253E1AAD6148D2uL), (+1, -1, 0x844B3BD5656E9920uL, 0x1BB5CAF854AB6F65uL)),
                ((+1, -6, 0x8053A67E02C366D8uL, 0xA517E790EDC29A55uL), (+1, -3, 0x956B553B4DF2529DuL, 0x5C39E4B0069BA00AuL)),
                ((+1, -9, 0x895FF2A3E256FBA1uL, 0x5695AAA42714E8CCuL), (+1, -6, 0xE1A3B33EF93AFD03uL, 0x231E6C5D110699B2uL)),
                ((+1, -13, 0xB05E1548E7E23D48uL, 0xFB5AD4FF1F622158uL), (+1, -9, 0xE0BE879EC29770EEuL, 0x507DA7B2479EFCB4uL)),
                ((+1, -18, 0xDAC30585EFAC39F7uL, 0x9DA0E8DACF612ABCuL), (+1, -12, 0x89C6EEA2C2AD4F43uL, 0xEAD70D9B3F2AED87uL)),
                ((-1, -27, 0xD8120D16A9E76F2FuL, 0xD6440972A142FF10uL), (+1, -17, 0xA62E97710E40ED61uL, 0x1BA02F07CF8E9B16uL)),
                ((-1, -27, 0x8A026CA95037830AuL, 0x6CDB8616BF9B41B1uL), (-1, -26, 0xC616DF3AE54FF8D0uL, 0x8BAE5EF25D84C1A8uL)),
                ((-1, -33, 0x9236DE9266373944uL, 0xC05046E9136C2EA8uL), (-1, -27, 0xD557AC4B35669014uL, 0x17475B15E0C25054uL)),
                ((+1, -49, 0x97FEF64008D2C4F5uL, 0x4A5B8122A53F792BuL), (-1, -33, 0xE4899279B46EA85BuL, 0x138432F7F5270CB2uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xB150AD2AB308BBC3uL, 0x53843EC3E586C405uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xB86BE0C61270AE62uL, 0x56B4CD2736AF54C8uL), (+1, 0, 0xA17D2B829E4CAB41uL, 0x7EB7BEFD42A9316DuL)),
                ((+1, -4, 0xAB1ED9F824C7CF10uL, 0x547F17F9A562B423uL), (+1, -1, 0xB7E703E51CC179B9uL, 0x9B97868EE8A12F92uL)),
                ((+1, -6, 0xB95023CCB96F39B9uL, 0xAA69F09C0E51622EuL), (+1, -3, 0xF8561E3D7D413DFBuL, 0x3FB9A8B995297B84uL)),
                ((+1, -8, 0x80A91070D8F828C1uL, 0x5728EA22D75CCB82uL), (+1, -5, 0xDBA7EB76D8675E82uL, 0xAD038DEB564B9EA9uL)),
                ((+1, -12, 0xECED065C0A414D78uL, 0x9A46167EBAA97BDAuL), (+1, -7, 0x84A9EB367531BD26uL, 0xDD7377766EC39963uL)),
                ((+1, -15, 0x908300F8B2FE7283uL, 0xABE107C3276D61CFuL), (+1, -11, 0xDD353AE1130A36FBuL, 0xB0C42C2A710288CEuL)),
                ((+1, -20, 0xE25F66F25C5B00BBuL, 0x872FE0F4118F975CuL), (+1, -15, 0xFB59AB9A74A6C35BuL, 0xB780BE8E14C2D288uL)),
                ((+1, -25, 0xD453399B910387B0uL, 0x2CB33FFB2A3E994EuL), (+1, -19, 0xBB5EF11E7601E0B4uL, 0x127D5323A0A46734uL)),
                ((+1, -31, 0xCFBE4C06412D7158uL, 0x5ACB93F9C2BEFE35uL), (+1, -24, 0xAA36A574D90780E8uL, 0xEEE95074702E3317uL)),
                ((+1, -38, 0x990D40CBBBF3065BuL, 0x201A8406514248F4uL), (+1, -30, 0xA3A5B2A75A3FD9CAuL, 0xC3F7409F84D1CDB0uL)),
                ((-1, -55, 0x85C514DD05D14AFCuL, 0x1C8F33E708543981uL), (+1, -38, 0xEFC9BE3E84E199F0uL, 0x7DE2193FE2878ACCuL)),
                ((+1, -64, 0xA444DCB5C5DA65EAuL, 0x95D03C3C327ED86DuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xB65E3A94FBFFA4E4uL, 0xE83D88483FA490F5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xF9004DB20DC19B71uL, 0x5BB6452496148EB3uL), (+1, -1, 0xCF27DBFA77DAD193uL, 0x377C8A9BB5AADD91uL)),
                ((+1, -6, 0x9358703BC93B550AuL, 0xF9062D5E1DC142B4uL), (+1, -2, 0x9414B7E4E734AC30uL, 0xF8AD0FD52AC5F554uL)),
                ((+1, -9, 0xC67FFFD4A4EC9F43uL, 0xAF4B8B0C25D0DF32uL), (+1, -5, 0xF6202F2E6F80D43CuL, 0x24E114014B497DF5uL)),
                ((+1, -12, 0xA7B40A0D9AC69D91uL, 0x3998986B457CA8A1uL), (+1, -7, 0x838B3440470E127FuL, 0x3A65BCE2F215DDB3uL)),
                ((+1, -16, 0xB86177326D668893uL, 0x0AB9547115971FDAuL), (+1, -11, 0xBCD02D265EB4A6F8uL, 0x926CA058049C52B5uL)),
                ((+1, -20, 0x8466CB4B6D117018uL, 0x65CAD742A26089AAuL), (+1, -15, 0xB85444F85305B48BuL, 0x1FE9E6F4959DC461uL)),
                ((+1, -26, 0xF2DFDEE5FCD53AE9uL, 0xDAA3D4BD3CD96F61uL), (+1, -20, 0xF2A2339C54D4AF85uL, 0x6930884CFD9EBD8CuL)),
                ((+1, -31, 0x87279FCFC53BA73FuL, 0x56F83AC6FAF15D76uL), (+1, -25, 0xD0F980483E076342uL, 0x45482BB778FAD9C7uL)),
                ((+1, -38, 0xA644EB02AB120C93uL, 0xFBE290CB35B327AFuL), (+1, -31, 0xDEAD230404A81B52uL, 0x6513B042937B3D00uL)),
                ((+1, -46, 0xBDF8606CC65EAD41uL, 0xE1AA483040BAB28EuL), (+1, -37, 0x8538312AE20B7B2AuL, 0xB2C21E28639E26D3uL)),
                ((+1, -55, 0x8724CA03222E7271uL, 0x0C6C05D833ACE7B5uL), (+1, -45, 0x95EEA0A52B0D1E53uL, 0xF1CFE743EE3509B8uL)),
                (Zero, (+1, -55, 0xD431FA0F5CF374C3uL, 0xCF7202881882CF69uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xB23D8BE264FF9D8BuL, 0xA01DF12D0794F0C1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -7, 0xEA471A197F0E15E3uL, 0xEB3CB3645802EA87uL), (+1, -2, 0xC9A3A6141E92A2DEuL, 0x570D655B3A9BC11FuL)),
                ((+1, -9, 0x851A0ED25B5D5561uL, 0xC39826C338478508uL), (+1, -4, 0x8BB2647928020A9CuL, 0x7542D2157279F502uL)),
                ((+1, -13, 0xAAF89199BB751047uL, 0x8714A4BC97998E1EuL), (+1, -8, 0xDF85657A970364B5uL, 0x34C9FE25CB63C085uL)),
                ((+1, -17, 0x884A872A0B401677uL, 0x3E1E6BD6DE573566uL), (+1, -12, 0xE3EB360ECE5D2970uL, 0x5004803252CDC7C0uL)),
                ((+1, -22, 0x8B6AC2AD97D587ABuL, 0x63A70C48AD55A4F7uL), (+1, -16, 0x9A3D079446324D8EuL, 0xCFC719CD9B226041uL)),
                ((+1, -28, 0xB6E59DFBA7FD78ABuL, 0x66D90F64FC00EF83uL), (+1, -21, 0x8BED7DEA5068CB8FuL, 0x6E56A42AF7CF9AFEuL)),
                ((+1, -34, 0x95617E907BB1205CuL, 0xC3451C369DFC90C3uL), (+1, -27, 0xA7FCEAB582F9E37DuL, 0x795FE18F6A545C67uL)),
                ((+1, -41, 0x8E393B8957BDDCBBuL, 0xA7698C5D514759DBuL), (+1, -33, 0x809B1DC7CF8CB2BEuL, 0xA9229926C3D46E04uL)),
                ((+1, -49, 0x8A556708BE148A8AuL, 0x3BBBCF682CF06AA7uL), (+1, -41, 0xE9F8D2119A7DE3AAuL, 0x700BDB1C2838D4BFuL)),
                ((+1, -59, 0xC7F71FC5FC58717EuL, 0xAAE84886E61F58B4uL), (+1, -49, 0xDCDBA7F1FC8537DFuL, 0xB90C6624CD701225uL)),
                ((-1, -80, 0xA5E7D96984317705uL, 0x1306C0D0B8ABB2E8uL), (+1, -58, 0x9CF3D8CA1FF7426AuL, 0xCA57AFE6F99742BDuL)),
                ((+1, -91, 0xC6171D0821DC9E79uL, 0xEA84D222525FEB51uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0xAD0EF8CB1CC868AFuL, 0xC54408D4D3BCBC00uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -9, 0xDFD686F4589B93D0uL, 0xF8EBBFD368AA2B18uL), (+1, -3, 0xC6CC9D7039241585uL, 0x7987CCF25347FE47uL)),
                ((+1, -13, 0xF95EC3CCE0E8339EuL, 0x0B0899BF5317E539uL), (+1, -6, 0x8750778E587A1E9FuL, 0xC3CFE14735D72136uL)),
                ((+1, -17, 0x9C7EE1F0C80596C2uL, 0xE3E5058C90FAF3FDuL), (+1, -11, 0xD3FCB1FBA12AFBDCuL, 0x023DB3EB5733EF2BuL)),
                ((+1, -23, 0xF310A797D9CCE061uL, 0x8EFBD9C6FE100CD0uL), (+1, -16, 0xD2F585527ED018C2uL, 0x862678E1A00D4F0BuL)),
                ((+1, -29, 0xF1A090542A678185uL, 0xB1BEDA5575B18ECAuL), (+1, -21, 0x8AEAD42B133181E6uL, 0x4700136EE68F080DuL)),
                ((+1, -35, 0x99BDBE94E2CB1E87uL, 0x0BDFF6F9F3AB9E9EuL), (+1, -28, 0xF4A6014BDEB8FED8uL, 0x0B121C9F102FD8AEuL)),
                ((+1, -43, 0xF355DBD9A624CCCAuL, 0x0CB4B80CAE70ED5DuL), (+1, -34, 0x8E414B6D618A1471uL, 0x576AA42A4BFE4F9BuL)),
                ((+1, -51, 0xE063735153F28CF1uL, 0x2470AD41ACEFE719uL), (+1, -42, 0xD2B587946091BAE3uL, 0xC5F12D5774E476A6uL)),
                ((+1, -60, 0xD371B947AA087E9AuL, 0xEE453867C40977E9uL), (+1, -50, 0xB94FD41C9DD487D1uL, 0xC04A1DF30C8F5A3DuL)),
                ((+1, -70, 0x9432B21781A14671uL, 0x4454A556EA704830uL), (+1, -59, 0xA927AA14304B57B9uL, 0x455EA09BF978CDC5uL)),
                ((-1, -94, 0xE84AB38137FE2023uL, 0xCB96DE744543B8CFuL), (+1, -70, 0xE8B78D3CDD023F0EuL, 0xC02A41831E64BA6AuL)),
                ((+1, -105, 0x870DF1C12B177D18uL, 0xA6D2254705617223uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA927397081BB773FuL, 0xF18FBF39DDBBC790uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x9E850E4D7E20CCE4uL, 0x9D9F394C7F3B988CuL), (+1, -1, 0xF4B3FA883C8895A3uL, 0xEC9BC17ACC7779F3uL)),
                ((+1, -2, 0x8B35AD75617C2960uL, 0xA79824A28471DFC9uL), (+1, -2, 0xD9BDE09EBA1335E5uL, 0xD9766DE56D39A4E8uL)),
                ((+1, -4, 0x95CD639831C409F4uL, 0x4ED7443601E166D8uL), (+1, -4, 0xEBAC4217ADD7B53FuL, 0x6FE063445498C89EuL)),
                ((+1, -7, 0xDC0E42D53640F444uL, 0x8F524CCBA4FDB0ABuL), (+1, -6, 0xAD191D4157EF8C8EuL, 0x6DB82991E394B2E0uL)),
                ((+1, -10, 0xE9F18C47417D39E2uL, 0x954D48A80D2737ABuL), (+1, -9, 0xB7988B491D363DCDuL, 0x3F16920C04505A24uL)),
                ((+1, -13, 0xBA76E8FD0E444DAFuL, 0xEBF362F73253618EuL), (+1, -12, 0x9256400DB6BBA06DuL, 0x8C14C2F43C2BCB73uL)),
                ((+1, -17, 0xE2FBF1E1C36D0FBAuL, 0x56A389AED0EA438FuL), (+1, -16, 0xB275252185B6FAB4uL, 0xC83FEDA6C8DB75BAuL)),
                ((+1, -21, 0xD33D30F9F6947331uL, 0x3C40D4E25031CB0CuL), (+1, -20, 0xA5E8D42FA75B7693uL, 0xC5D27641084E1B43uL)),
                ((+1, -25, 0x93554E98A1E46791uL, 0x830B9C652DAE004EuL), (+1, -25, 0xE717C46A76235571uL, 0x9FF090A5D5C3FEBDuL)),
                ((+1, -30, 0x8FCF3C1B4DC6C2D6uL, 0xD5924BD292C81806uL), (+1, -30, 0xE26CECB05DF09C94uL, 0x9D55A6A1493CDE6AuL)),
                ((+1, -36, 0xAA376EDACF4E5056uL, 0x2F708901D7E51707uL), (+1, -35, 0x85684CA45B8326A2uL, 0xF1DACD59A83E6245uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA511C4AA57103B34uL, 0x483DAA3CDB386288uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xF2862C6CA1F58670uL, 0x5A674C5440860B92uL), (+1, -1, 0xBDE3146A3F2BC1E9uL, 0x34653A9A837E8365uL)),
                ((+1, -3, 0xAF3CFABE62C7A8D8uL, 0xEA1BC863DB5BA3B1uL), (+1, -2, 0x899F3E7DF0FA3CD9uL, 0x2C7F8F1CA0FB8C14uL)),
                ((+1, -5, 0xA38D82C0CA1A58CAuL, 0xAC1D979F382DD51BuL), (+1, -4, 0x807E0CDA72DE719AuL, 0x1317585809EED6CCuL)),
                ((+1, -8, 0xDA7393141670036CuL, 0xA9AAF7007EDF0C18uL), (+1, -7, 0xAB99EAB7551512F8uL, 0xE56975C7DC9C045AuL)),
                ((+1, -11, 0xDB206FC9546DBB55uL, 0xA0F622CAF402F6A9uL), (+1, -10, 0xAC1549E816CE2954uL, 0x01388EF77B2C749FuL)),
                ((+1, -14, 0xAA83EC15177DDDB9uL, 0xC7E54C689C18EAE6uL), (+1, -13, 0x85E482A9553B1872uL, 0x8F4E008B6C525D22uL)),
                ((+1, -18, 0xD405FA9BB65EA116uL, 0x5FEC79B23B3AA985uL), (+1, -17, 0xA692120982F4F6D0uL, 0xFB2088BAE760D743uL)),
                ((+1, -22, 0xD93A3066CEAB79DAuL, 0x68B4CE627918BF90uL), (+1, -21, 0xAAA5764F6DF04608uL, 0xCB431F96DEE5167AuL)),
                ((+1, -26, 0xBD855A5279F944AAuL, 0x793274B82C0D9C69uL), (+1, -25, 0x94AB6DF9DAB694F5uL, 0xCBFB223D127EC53AuL)),
                ((+1, -30, 0x920A0ACB713BF054uL, 0xA8C9F5D646B57244uL), (+1, -30, 0xE5E9689E44670CADuL, 0x286BCA24DD6D3640uL)),
                ((+1, -35, 0xCBB7E376125D65B7uL, 0x7638EF26AC8DD5BFuL), (+1, -34, 0x9F8FF76989870143uL, 0x1D05B8171D333738uL)),
                ((+1, -39, 0x81DF29070D38167CuL, 0x120812329873A793uL), (+1, -39, 0xCC80975F73E93EF5uL, 0xB2630EE6F6EFE59EuL)),
                ((+1, -44, 0x918690F53216B956uL, 0x7459C80B6057D99EuL), (+1, -44, 0xE432BD19E900F57CuL, 0xC5050A1B8F97B9C2uL)),
                ((+1, -49, 0x8E9A62D7289DC528uL, 0xD0B3573092865D26uL), (+1, -49, 0xE034E4A0C7226358uL, 0x526A21C1D954F8F2uL)),
                ((+1, -55, 0xC6F647FFF8D22F71uL, 0xA05A4650FE136D2FuL), (+1, -54, 0x9C32D50AC0E429FCuL, 0x62BC6AFF879403D9uL)),
                ((+1, -61, 0xFA6F76B1EDE2D7A0uL, 0x0A74FA6A1DA0A35BuL), (+1, -60, 0xC4B6064588D99B93uL, 0x069DFB7A059DB6CFuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2FDD652CE522608uL, 0xCFE69252DDD1A0B3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x86BF4ECCACB51B21uL, 0xF053CE459C094876uL), (+1, -2, 0xD3AC0547D9576AECuL, 0xCE9F5830957A3E60uL)),
                ((+1, -5, 0xD7747CE10A5B04C9uL, 0x0246373FE25E1C2EuL), (+1, -4, 0xA936D1884129BA65uL, 0x439C50C844B832C4uL)),
                ((+1, -8, 0xDD982CBBB043E959uL, 0x5C42DA5037C88832uL), (+1, -7, 0xAE0A6A89906E272FuL, 0x2FB4C535B2F3C9B3uL)),
                ((+1, -11, 0xA45EA0A544D4987EuL, 0x75D7DD728A622D56uL), (+1, -10, 0x81189E70AADE3F06uL, 0xB1C62F8798662CA5uL)),
                ((+1, -15, 0xBAC03A40FF73A735uL, 0x12F0888699DEA336uL), (+1, -14, 0x92AC3635E67C28DCuL, 0xBFC7BB84C9F63506uL)),
                ((+1, -19, 0xA84DBBD36B3933D6uL, 0x3A6E78667832E74CuL), (+1, -18, 0x842FB31860260A0BuL, 0xDD30AC4276D69BE1uL)),
                ((+1, -24, 0xF5A95A9F4677474BuL, 0x80AA869F25FA5672uL), (+1, -23, 0xC0F0EFE6A06D2FBCuL, 0xEA5EFDE5761774A2uL)),
                ((+1, -28, 0x92C2F5EBE5DF84B1uL, 0xEC4F15FCA9CCD0ABuL), (+1, -28, 0xE688723049A30295uL, 0xCBAC45FAF25A268DuL)),
                ((+1, -33, 0x8FDFE106555E2F95uL, 0xBB5C7281E8DB3557uL), (+1, -33, 0xE1FF6F990B36F6B2uL, 0x0E99DAC509B19617uL)),
                ((+1, -39, 0xE5ECFC03C2BBFE4EuL, 0x897B2843ED992A8DuL), (+1, -38, 0xB49524F9486C32E1uL, 0xAF9ABD88DCFF2A0CuL)),
                ((+1, -44, 0x92F6E9E0AB231AA8uL, 0x30CFE9861F08644AuL), (+1, -44, 0xE6DA44C1E117A310uL, 0x5A1E44280399429BuL)),
                ((+1, -50, 0x9097DD3644820B4AuL, 0x3A5CCBE8AB36437AuL), (+1, -50, 0xE31FF45C2295197EuL, 0x5DA606296B33F761uL)),
                ((+1, -57, 0xC92070967BC5DB68uL, 0x4A986E9A69B01A41uL), (+1, -56, 0x9DF74D8354A06C01uL, 0x40334E09D51D3E40uL)),
                ((+1, -64, 0x9FD69E276A95296AuL, 0xE30BDFA778D97552uL), (+1, -64, 0xFB11EC44AFD18594uL, 0xFFA586098B7CF7ACuL)),
                ((-1, -88, 0xD353FEC5B6EC84DDuL, 0x38A0F6324C09D723uL), Zero),
                ((+1, -96, 0x90010FEFD8879238uL, 0xDA581993B35AE0EBuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F983771FB6B168uL, 0x8E3AC18FE654D98FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF4E2F45962FECF02uL, 0x4E6665695BF0D36FuL), (+1, -2, 0xC0555C52C55D9E1BuL, 0xE437C05898F28AA4uL)),
                ((+1, -5, 0xB4BD1A15754518CBuL, 0xE888AB54FC817678uL), (+1, -4, 0x8DF3AAE68C3C3747uL, 0x0AB78EBD599B2C69uL)),
                ((+1, -8, 0xAE31199CD200568FuL, 0xBFD42820D5CDDDBAuL), (+1, -7, 0x88CF56C3C411E1FAuL, 0x9D68C8BC53C0EE16uL)),
                ((+1, -12, 0xF5D57C9E5A28B93FuL, 0x43ADB89BBF278F97uL), (+1, -11, 0xC113D84B274162DCuL, 0x768CF78DD627BD7FuL)),
                ((+1, -15, 0x86FD0C806F9A6AFCuL, 0x5F1DC799FDA43AC9uL), (+1, -15, 0xD40A15DD515D14F0uL, 0x870F14C002FE3DF7uL)),
                ((+1, -20, 0xEF55EB9EB706B5DCuL, 0x6626FBC1DA072F95uL), (+1, -19, 0xBBF9486DE12D2A33uL, 0xABB5225CA3497ABDuL)),
                ((+1, -24, 0xAF4C8D3EA61BCF70uL, 0x7D282A87E4D68E6DuL), (+1, -23, 0x89ADF61944DD282CuL, 0x4CCB24FAA917904AuL)),
                ((+1, -29, 0xD749F283AF63F495uL, 0xCD44C5AE431BA6C7uL), (+1, -28, 0xA916649B7AA89A1BuL, 0x333F01AA8E95972AuL)),
                ((+1, -34, 0xDF7FB39B0CBF059DuL, 0x5A957EC5640E7785uL), (+1, -33, 0xAF891B60D1EC35A1uL, 0x32BB69F5C78DF932uL)),
                ((+1, -39, 0xC4A7F528729817F2uL, 0x9DBC6C9EF5B7E202uL), (+1, -38, 0x9A740D47914D7DDAuL, 0x1C3FDA4883C54CB3uL)),
                ((+1, -44, 0x9237D06F1A4B964DuL, 0xD8977356FDED497FuL), (+1, -44, 0xE5ADC199413F3235uL, 0x762CFEA311ABD687uL)),
                ((+1, -50, 0xB5E0F63A33CE36A0uL, 0x041A9FC041F71E69uL), (+1, -49, 0x8ED8E4D6DA562160uL, 0x1631EDE83F81459DuL)),
                ((+1, -56, 0xB95B839165AE6E8DuL, 0x848BC360885BA1BDuL), (+1, -55, 0x919454F733E2D675uL, 0x368317630E0A71C8uL)),
                ((+1, -62, 0x94BD99040D950D12uL, 0x2470C8A6659BF3BDuL), (+1, -62, 0xE9A4268F9904A4ECuL, 0x145773E3C9E2C004uL)),
                ((+1, -69, 0xACFE2AD9E0382B77uL, 0x1BB45E3160FD09CEuL), (+1, -68, 0x87DE464794297BB7uL, 0xCC507889734A8265uL)),
                ((+1, -77, 0xEAF68C5FD67E407FuL, 0x73045CCA3ACF44C4uL), (+1, -76, 0xB88A21543B13291AuL, 0xAB2D41BD0A843EBFuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E44153BuL, 0xCAEA066EAE6DA3D0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8EB071CC23A5F6CDuL, 0xBD1D171EA28E5952uL), (+1, -2, 0xE022BF192444C354uL, 0x839EEDBA7CB2F087uL)),
                ((+1, -5, 0xF60387AE85F3FFB4uL, 0x7B094E9B20BCDBDDuL), (+1, -4, 0xC13801D7CA48081BuL, 0x2C3E331AC54ECE74uL)),
                ((+1, -7, 0x8AFCB37BF1AEDA29uL, 0x57A66232DBF72C23uL), (+1, -7, 0xDA5208DC17D40C09uL, 0xAE9F1AB6666A7DA8uL)),
                ((+1, -11, 0xE72288A0DA4AFF38uL, 0xD8B0F34C4C7E7A3DuL), (+1, -10, 0xB5886DB2172F2571uL, 0x4807483F003A5B86uL)),
                ((+1, -14, 0x9686B785C870553DuL, 0x5A6FE0EAF9941B89uL), (+1, -14, 0xEC723115B0C7EEDDuL, 0xF845BAAB1C4E1A73uL)),
                ((+1, -18, 0x9F88CEAA2AB0A284uL, 0x8BA973B84B170A78uL), (+1, -18, 0xFA98970E81BD9676uL, 0x964E54FCEC9148ECuL)),
                ((+1, -22, 0x8D1004B4B7291B99uL, 0x0031F0FCED73ABD8uL), (+1, -22, 0xDD94A03649E75512uL, 0x4A191A0F2FA18D92uL)),
                ((+1, -27, 0xD391E33868980BEEuL, 0x7198D56EA56289E7uL), (+1, -26, 0xA62AA5A3D3D9B719uL, 0x4A84178BF6F1119FuL)),
                ((+1, -31, 0x8827793DDECB70DDuL, 0xF5844D1D3B240C23uL), (+1, -31, 0xD5DED9930D1800BCuL, 0xA9E40CC9AD7B5ADAuL)),
                ((+1, -36, 0x966E4F5BD6D5BA9FuL, 0xBB4C7EF2BFCDEF28uL), (+1, -36, 0xEC4BDA7DEE27348DuL, 0xF64182CEEF499D2DuL)),
                ((+1, -41, 0x923245C9CB03A4C8uL, 0x0E0AEC27AEA93B18uL), (+1, -41, 0xE5A50D31DB1620FDuL, 0xCC519494C1F00CD5uL)),
                ((+1, -47, 0xE3E2905BE3427CACuL, 0xDCC5529D78C37014uL), (+1, -46, 0xB2FB003DDCDC93FCuL, 0x9DDB0AC81B3EBF3EuL)),
                ((+1, -52, 0xB9600E667AC78953uL, 0x7D1863F4B5E77092uL), (+1, -51, 0x9197E64081A57B6DuL, 0x46B3F1BEA136C8C1uL)),
                ((+1, -58, 0x8758B5069B2C61B0uL, 0x7857D3EC3E8B0B74uL), (+1, -58, 0xD49A0FD81A59E157uL, 0x920C34729324A358uL)),
                ((+1, -63, 0x87FFF4297077967CuL, 0xC0D0B15288BC670CuL), (+1, -63, 0xD5A0C5B3DEBDDB92uL, 0x9CD42FA8870F60EAuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_0 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xC4F009B69250EDE7uL, 0xC9099B4F040486C3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x8D9FFE34F8D57E10uL, 0xDA8C4B7C0EDAC188uL), (+1, 1, 0xDC9A8AB8CCA6AE46uL, 0x3C0BA21F6DD94066uL)),
                ((+1, 0, 0xBAA19411DC1EBCA4uL, 0x30566ED17682CB44uL), (+1, 2, 0xC12FA0B1113C1A41uL, 0x5592AD9764CA3724uL)),
                ((+1, 1, 0x955FA72F33C1BD53uL, 0x032B9D0DED09250EuL), (+1, 2, 0xDE4FDB8477E67BDAuL, 0x3716C479ECB905CFuL)),
                ((+1, 1, 0xA26E582E11CA3884uL, 0xCEAAE277BCA83DDDuL), (+1, 2, 0xB9B166C2810BEF7EuL, 0xF0EE3A69981A3678uL)),
                ((+1, 0, 0xFE6D5CAB41EBD774uL, 0x45871C18A24491AEuL), (+1, 1, 0xEC448DE9C1B1AF8AuL, 0x97809768B930CB00uL)),
                ((+1, 0, 0x94758F2FE182BC5AuL, 0x27D18EF5BCCC613AuL), (+1, 0, 0xEA8D184F69A46DFAuL, 0x901C3BC7F33A9B82uL)),
                ((+1, -1, 0x8384E66515230239uL, 0xC040D6A2501CBFF0uL), (+1, -1, 0xB759310AABE21733uL, 0xA3EF22BFAB10D9C7uL)),
                ((+1, -3, 0xB1E20422AE7C8A94uL, 0x0B7B9238DE195634uL), (+1, -3, 0xE0FFCEB37A6A90D5uL, 0x69708DCCA49C260FuL)),
                ((+1, -5, 0xB5E50111AA856EDEuL, 0x0C4D12372FF2C403uL), (+1, -5, 0xD53BCCCC8FA334C6uL, 0x4A2CE76952180B35uL)),
                ((+1, -7, 0x88624179A025FAD8uL, 0xE71A76CADB093261uL), (+1, -7, 0x96A0FADA3086D470uL, 0x4FFBCD38096B3D07uL)),
                ((+1, -10, 0x8C62B9073228EF4FuL, 0xC4BBC44220140368uL), (+1, -10, 0x93F084B92159BD47uL, 0xE79AC272DBCE48D1uL)),
                ((+1, -14, 0xAA342AC4BE9FE6F4uL, 0xD17B79F02F3E4CFBuL), (+1, -14, 0xAD244B8C3509ACF7uL, 0x30535DA9D73ED6FCuL)),
                ((+1, -19, 0x92788AE2D61C06F7uL, 0xECA3F58BEFBD57BDuL), (+1, -19, 0x9277B79C88771461uL, 0x57D862EA8F03A883uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1p5_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0x9999B0DCC109C797uL, 0x824FE2A23D327F7CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x8B0B5B75A0F8EF0FuL, 0xE65434A1748289F7uL), (+1, 1, 0xAF3456DDC01954F5uL, 0x9929C9550CE65FB7uL)),
                ((+1, -2, 0xE456BA3AF2429857uL, 0x13DFFB1F7FEC19C1uL), (+1, 2, 0x8DF837A449227E6DuL, 0x5B8C7F50FB36912DuL)),
                ((+1, -1, 0xE0681C314A00B281uL, 0x754AA1492C010A63uL), (+1, 2, 0x9A655626D0786410uL, 0xA66A97D85830D845uL)),
                ((+1, 0, 0x92B6E148362EE85EuL, 0xE4807442F0440E3FuL), (+1, 1, 0xFD615DD15A0691E6uL, 0xCC1741BF5602C4BEuL)),
                ((+1, 0, 0x864C019028CC8379uL, 0x60DBD3F3E3A004DDuL), (+1, 1, 0x9F6CB3BAC32C0E45uL, 0x51A6CB3828349842uL)),
                ((+1, -1, 0xB0A570E2DC7359D4uL, 0x3A6B06F4D951043FuL), (+1, 0, 0x9F04C5E2B9314A98uL, 0xC0F5A3955C40D6ECuL)),
                ((+1, -2, 0xA9083A880BA1E2BEuL, 0xABA69022152F8104uL), (+1, -2, 0xF78FF5C01A3CC2F8uL, 0x54CF17E474217497uL)),
                ((+1, -4, 0xEBFD17F254EC5471uL, 0xDECCCEFCDE480818uL), (+1, -3, 0x96C19E94AC9EE834uL, 0x3416B7E39B21EF13uL)),
                ((+1, -6, 0xECE4EAFD421F2B7DuL, 0x6BC02303878ABB8BuL), (+1, -5, 0x871A1EF668F201CEuL, 0xDBAE1EDE933A460FuL)),
                ((+1, -8, 0x9F7308FF07B950F3uL, 0xE3728C3E76541410uL), (+1, -8, 0xA946572121818749uL, 0x85F1633671692561uL)),
                ((+1, -12, 0xDD3AEDD098C78D95uL, 0x8F5F6B0939DEE296uL), (+1, -12, 0xDBC7A200F8577613uL, 0x024F021EB70814E0uL)),
                ((-1, -23, 0x895294314A3A1B60uL, 0x4187B4673081C3B7uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1p75_1p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -8, 0x9F7562F6138D1570uL, 0xE743BB1700BBF576uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xABD9438AE4CABA33uL, 0xEB48754C59D6E173uL), (+1, 1, 0x8DB2AD9504BF96D2uL, 0xB596238B6142D927uL)),
                ((+1, -3, 0xA7F9C2B3EE41658CuL, 0xEB55DD61AA4EA3F9uL), (+1, 1, 0xE4DE79BF47D51060uL, 0xD26B5AB3A511A2F2uL)),
                ((+1, -2, 0xC4427A55A3D006A0uL, 0xFB769B273DE45950uL), (+1, 1, 0xEDC3DE19B9A6FF7AuL, 0x532091313BC2F496uL)),
                ((+1, -1, 0x980C246ADCA3531FuL, 0xB7D558092273951EuL), (+1, 1, 0xC9BBC68EE0B2F4A8uL, 0x9047527A8DB8429BuL)),
                ((+1, -1, 0xA3EAC8921289370DuL, 0xC9AD2D71E13599ABuL), (+1, 1, 0x801792B9D2404FC9uL, 0x1F70CBBB498C5DD9uL)),
                ((+1, -2, 0xFB7411957B014AE2uL, 0x22E7C224AA93FA47uL), (+1, 0, 0x89113AB5DB7A7BC5uL, 0x9968AA21F392E7A2uL)),
                ((+1, -2, 0x8A8883D511FEB2F8uL, 0x700A5323A819DBE8uL), (+1, -2, 0xDFB8AEE288564B01uL, 0x0716E6DCF2FE79C9uL)),
                ((+1, -4, 0xDCAE866E3AF8B90FuL, 0x25926C03FAE240A0uL), (+1, -3, 0x993558A07B432851uL, 0x8313E2C5ED502BF9uL)),
                ((+1, -5, 0x8006270AB193C84AuL, 0xAC1D38CDEAC2471FuL), (+1, -5, 0x969FDEBADF8825FDuL, 0xC61988FC551FF4F0uL)),
                ((+1, -8, 0xD2D0E374B03AE749uL, 0xE70B9617D7F6E052uL), (+1, -8, 0xE5D5702EEEC123F0uL, 0x83F9E404B4586881uL)),
                ((+1, -11, 0xC2FE9ADC3E9F0179uL, 0x0C6AEE155A0A1B0CuL), (+1, -11, 0xC127E2167BEBE843uL, 0x5FD28850BDB966D1uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_1p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -11, 0xB95D9EFD18468B3DuL, 0x3FDF7F5A8B822763uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -8, 0xFF640C1A9E8AF345uL, 0x844AA55BB14EEFA8uL), (+1, 0, 0xE8D675F892BA661AuL, 0xADDEB985BB853881uL)),
                ((+1, -5, 0xA19D6EBBA36105B0uL, 0x0FE3CB11AB0D24ABuL), (+1, 1, 0xDCF0C08689A0B11BuL, 0x1184F51026F75EE5uL)),
                ((+1, -4, 0xF81C1043FC28229CuL, 0xEB3EB6A4611D1A68uL), (+1, 1, 0xDBB3E547B44A8A00uL, 0x0CDC039E5C6E1767uL)),
                ((+1, -2, 0x8062C925D5456BB3uL, 0xD2823FDBB4252CEEuL), (+1, 1, 0xDF3250DEC117A98AuL, 0x5831E514219423BFuL)),
                ((+1, -2, 0xBC2A14E9915F440FuL, 0x20F53D126106DF37uL), (+1, 1, 0x8C0CC57F7F4E9383uL, 0xCC4BE45CFB9E911CuL)),
                ((+1, -2, 0xC7504C065EAEF3ACuL, 0x03ABA0BDC60D5E0AuL), (+1, 0, 0xBBCC1FA7C67C94C1uL, 0x3D0044583D00F808uL)),
                ((+1, -2, 0x98940288B6F54710uL, 0xD45A8C07F01E554FuL), (+1, -1, 0x92DC7DB264728F1DuL, 0x1BB32B736AFDCAFBuL)),
                ((+1, -3, 0xA540374E5A91AD03uL, 0xDCAFFCCE858C0F47uL), (+1, -2, 0x8A6F5586D173CA11uL, 0x1CBCDDA630702DC3uL)),
                ((+1, -5, 0xEFE21E1725F4B842uL, 0x3671CF32FBD4FA15uL), (+1, -5, 0xD9C4DA8874B5F534uL, 0xF0BD7905AD0F9ABBuL)),
                ((+1, -7, 0xCD156F0B88DEA391uL, 0xAC8F9F91F74E1AEDuL), (+1, -6, 0x93B71EA3D9DC1C2BuL, 0xAB6D4D0A88179095uL)),
                ((+1, -10, 0x82D1E8D0257CB937uL, 0x0AD74842D3172665uL), Zero),
                ((-1, -14, 0xA82CA95369EECFBDuL, 0xA8C67E226399D771uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xBE4A10B98B82AAF6uL, 0xD2A6D823E04C1F29uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x8E0BA20FA39923EEuL, 0xB59E87003E89CED9uL), (+1, -2, 0x8E4EE8F1D82E78F1uL, 0x8B6F76C73F43EB55uL)),
                ((+1, -3, 0xE465DBF7D9610E68uL, 0x2FE70E896B9C1797uL), (+1, -1, 0x9FF1EDFB87B7E4EFuL, 0xDFFC2E4440043BDAuL)),
                ((+1, -4, 0xB01A113D8C3E4769uL, 0x6B1C889CEE37DC3DuL), (+1, -3, 0xB3C971A1F79050A8uL, 0x7CC7F0657733CED5uL)),
                ((+1, -5, 0xF1061D2B15FB7EDEuL, 0x36D5734636C25841uL), (+1, -3, 0xAC3A4A1CBC97D04BuL, 0x5E9C33611C120CF6uL)),
                ((+1, -6, 0xAEC29ADED6EE3754uL, 0xFDD82D2543E1A4B0uL), (+1, -5, 0xB8152E856E7A6B07uL, 0x5DEE75CE878956B3uL)),
                ((+1, -7, 0x8F7A2D261C74A744uL, 0x598A379F1FF22F85uL), (+1, -6, 0xCBFBBCA2A7EA89C4uL, 0x7C9A4AFBA27F8E60uL)),
                ((+1, -9, 0xB1CA0B53C0A33A58uL, 0x46EB6270D1B89174uL), (+1, -8, 0xC2E1CD099C99A069uL, 0x0BE9BB51F1943C8BuL)),
                ((+1, -11, 0xC9AA83AF3AFD073BuL, 0x95D703E4F0FEDAE3uL), (+1, -9, 0x8BC91209F1D08165uL, 0x8992BDC68C157500uL)),
                ((+1, -13, 0xC2D5ED7144D19599uL, 0x8EF3D05D227825C0uL), (+1, -12, 0xDEE1A143AF694FBAuL, 0x6EF7F48428A6BD50uL)),
                ((+1, -15, 0x9E8CD3842B0685BFuL, 0x70ECBD860203DFE6uL), (+1, -14, 0xD46E4D35A9489BDAuL, 0xD7531F48E94C6E5BuL)),
                ((+1, -18, 0xD86DF563E54B194FuL, 0x9033B99B8B542B00uL), (+1, -16, 0x80C110EE6E0AEEB4uL, 0x67982BA86E0856C4uL)),
                ((+1, -21, 0xECA1F634323BE600uL, 0x2D9490F281A0D057uL), (+1, -19, 0x9962BB8DFE339285uL, 0x02ACD1B3763E8BD7uL)),
                ((+1, -24, 0xC1DC47806A64F434uL, 0x3FC65D963AFAA37CuL), (+1, -23, 0xEDB9945EF64908C9uL, 0x2474A2875074A9BDuL)),
                ((+1, -28, 0xD55A18C75414DB29uL, 0x62731A45558F0470uL), (+1, -26, 0x875E4FB567E4CA90uL, 0xEBEE4F8B18987F8AuL)),
                ((+1, -33, 0xED46C91F07461998uL, 0x237479639F523E9EuL), (+1, -31, 0x93B1DC8E5831C57CuL, 0x329DBA3161770E7CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xCB853D0C0B2B6925uL, 0x3BEC93BF300864A7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xD2E29085E3221CBCuL, 0xD0255D951C241C9DuL), (+1, 1, 0x84457891E4F26B0AuL, 0xF0457E76CEB09B88uL)),
                ((+1, -1, 0xE03683E92F5C1CAEuL, 0xB21CCCB070BCF0EAuL), (+1, 1, 0x8C8BC34715B06D57uL, 0x6410105E3F19EAC9uL)),
                ((+1, -1, 0x9C421E920FB40C7DuL, 0x39EEA3680FEAA760uL), (+1, 0, 0xC3D5B85967E95DDCuL, 0x9B23991D75FA22BDuL)),
                ((+1, -2, 0x9D8EB29B1FD3766AuL, 0x75DA614F4D88AF30uL), (+1, -1, 0xC576B15C90D17393uL, 0xA948FC0BFBFCD98DuL)),
                ((+1, -4, 0xF160904744C0C47DuL, 0x55576F5A50D11AEFuL), (+1, -2, 0x9741D29D7A0043B5uL, 0x845A96D51467D723uL)),
                ((+1, -5, 0x909A398CEEE03AB9uL, 0x8C8C0F1E61057AB0uL), (+1, -4, 0xB53CC67317B55409uL, 0xE7D04E317D414D32uL)),
                ((+1, -7, 0x89FD59D811774BAEuL, 0x0CC9654C0772BAE9uL), (+1, -6, 0xACF15C3233C1549CuL, 0x563B709ED6168D66uL)),
                ((+1, -10, 0xD46D42145066E0CFuL, 0x9E468D8B7735310FuL), (+1, -8, 0x851EA80FBB1E1453uL, 0x50B6BAA6ED6CEC25uL)),
                ((+1, -12, 0x850698728665E8E3uL, 0xD0EA75DF890258D3uL), (+1, -11, 0xA6B7C86EFF84A99EuL, 0xF1C2DF0CD55AD8CAuL)),
                ((+1, -15, 0x88391C366DDF10D2uL, 0xDC31C7CFFDCBD0A5uL), (+1, -14, 0xAABDDB346D0A3332uL, 0xAA651D85732C7055uL)),
                ((+1, -19, 0xE42B31B5C11A14E7uL, 0x28A3BB9136847DBCuL), (+1, -17, 0x8EF8361F7779C5C9uL, 0xFD7BF1C7EDDD4BFEuL)),
                ((+1, -22, 0x9B614FA7FAF560F1uL, 0xC902E17D3204420FuL), (+1, -21, 0xC2C340CF9745C8E4uL, 0x627495546DDEBDEEuL)),
                ((+1, -26, 0xA9868F80C33970F8uL, 0x6945BC1797DCEAF4uL), (+1, -25, 0xD4718D7F3A4E53F1uL, 0xFF0275ADE54E8446uL)),
                ((+1, -30, 0x8F0596EC2A4C0092uL, 0x8D8767D5E171AA14uL), (+1, -29, 0xB3456148CF06585AuL, 0x872CAD9C3C973B1DuL)),
                ((+1, -35, 0xAE2B8761CC54B00CuL, 0xD2AB02D49CB64BC2uL), (+1, -34, 0xDA44DF56A79BFF6AuL, 0x75D25E17CDA76247uL)),
                ((+1, -40, 0x80E4946F7F49F630uL, 0xA609B75C1B75761BuL), (+1, -39, 0xA18E7EA246169D30uL, 0x55593BDCDC089294uL)),
                ((+1, -61, 0xD159F532F718DEE5uL, 0x023A91A409B5673DuL), Zero),
            }));

            private static ddouble PlusValue(ddouble x) {
                Debug.Assert(x >= 0);

                if (x <= 64d) {
                    ddouble y;

                    if (x <= 1d) {
                        Debug.WriteLine("pade minimum segment passed");

                        y = ApproxUtil.Pade(x, pade_plus_0_1);
                    }
                    else if (x <= 2d) {
                        y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                    }
                    else if (x <= 4d) {
                        y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                    }
                    else if (x <= 8d) {
                        y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                    }
                    else if (x <= 16d) {
                        y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                    }
                    else if (x <= 32d) {
                        y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                    }

                    return y;
                }
                else {
                    int exponent = double.ILogB((double)x);

                    ddouble v;
                    if (exponent < 8) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -6)), pade_plus_expp6_8);
                    }
                    else if (exponent < 16) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -8)), pade_plus_expp8_16);
                    }
                    else if (exponent < 32) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -16)), pade_plus_expp16_32);
                    }
                    else if (exponent < 64) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -32)), pade_plus_expp32_64);
                    }
                    else if (exponent < 128) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -64)), pade_plus_expp64_128);
                    }
                    else {
                        v = 2 * RcpPi;
                    }

                    ddouble y = v / x;

                    return y;
                }
            }

            private static ddouble MinusValue(ddouble x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 2d) {
                    ddouble y;
                    if (x <= 1d) {
                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0);
                    }
                    else if (x <= 1.5d) {
                        y = ApproxUtil.Pade(1.5d - x, pade_minus_1p5_1);
                    }
                    else if (x <= 1.75d) {
                        y = ApproxUtil.Pade(1.75d - x, pade_minus_1p75_1p5);
                    }
                    else {
                        y = ApproxUtil.Pade(2d - x, pade_minus_2_1p75);
                    }

                    return y;
                }
                else if (x <= 8d) {
                    ddouble v;
                    if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }

                    ddouble sigma = Exp(x * pi_half - lambda_bias);

                    ddouble y = v / Sqrt(sigma) * Exp(-sigma);

                    return y;
                }
                else {
                    return 0d;
                }
            }

            public static ddouble Value(ddouble x, bool complementary) {
                if (x >= 0d) {
                    return complementary ? PlusValue(x) : 1d - PlusValue(x);
                }
                else if (x <= 0d) {
                    return complementary ? 1d - MinusValue(x) : MinusValue(x);
                }
                else {
                    return NaN;
                }
            }
        }

        private static class QuantilePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 2, 0xB5CFBD01F2E52D51uL, 0x63E7A633C6DD64B2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 7, 0xA85C7B426B738E09uL, 0x852B78B74132FBB4uL), (+1, 5, 0x97AD03F603B08CB8uL, 0x4F8637C6F287FF7BuL)),
                ((+1, 10, 0xF55990C458FE4BD2uL, 0xA211374EE80DBDEEuL), (+1, 9, 0x95A2106DCFB6EDC8uL, 0x9902FE027EA9FF3EuL)),
                ((+1, 13, 0xAAE8329FE2EDFDD5uL, 0x3E41D0DB8EA71AAEuL), (+1, 12, 0x9DA772E89EC9E41CuL, 0x41913E5DC4881C81uL)),
                ((+1, 14, 0xC214DBBA2F3CF708uL, 0x4256ABC678301B41uL), (+1, 14, 0xBB080F92AA21EB70uL, 0xC92F011A93090660uL)),
                ((-1, 14, 0x9EDADAB0E1700DE1uL, 0x518B46D55981A205uL), (+1, 15, 0xE78334802B4841C8uL, 0x18BECD671818FC86uL)),
                ((-1, 17, 0xBAEC2A4BA50F6589uL, 0x79791B69F34963FFuL), (+1, 15, 0xA7C51A2BDB85BFE6uL, 0x7DB24436E59EA244uL)),
                ((-1, 17, 0xB976AE9DE5F36F4FuL, 0x15FDDFEAACB0FE4AuL), (-1, 16, 0xDBEED9EDC85F4FB8uL, 0x4A62126838721F6FuL)),
                ((+1, 18, 0x9519CFF3E9025C3FuL, 0x3435CEF4F3046215uL), (-1, 17, 0xCA3B19B8B72EEEB7uL, 0x6CAA72148740723AuL)),
                ((+1, 18, 0xD141C526560DC83CuL, 0xAA73132EA59C544DuL), (+1, 14, 0xE4290C4306C230EAuL, 0xAE190AF9B1E74BE6uL)),
                ((-1, 17, 0xC5DBDEC1E3B9FA3DuL, 0x4BD9FB8248383F80uL), (+1, 17, 0xB498E3AFCBEBCB76uL, 0x7714D9A5C6DBF5B8uL)),
                ((-1, 17, 0xE164F376B61AF1DAuL, 0x9B0DB1360B659AF5uL), (+1, 13, 0x85711EFA71DE3311uL, 0x3DE81F11A76E4085uL)),
                ((+1, 15, 0xE3890B13CE8A2B74uL, 0x9A950CD04DB0C3B9uL), (-1, 15, 0xA14153D8F7B5D920uL, 0x0A596D2FBABD4C09uL)),
                ((+1, 13, 0xEFAE712565091A75uL, 0x982519362433AE68uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 1, 0xA340906BF59E3687uL, 0xF240AD51DB367076uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 4, 0x9B0A5F5E319BDC86uL, 0x110C1ECBFBE158F3uL), (+1, 3, 0xD0289F5B91A79EE5uL, 0x3FDF3F1BB1926220uL)),
                ((+1, 3, 0xF65E0BFF6BC927B7uL, 0xA8343138FEF6933CuL), (+1, 5, 0xE952BBB04B5704A1uL, 0xDFB5D2E0CEC6B301uL)),
                ((-1, 7, 0xB5FA4E2CFFACA121uL, 0x4E1021CF51829BB4uL), (+1, 6, 0x9CAB7EDB0CEF768BuL, 0x96473D2DA2BAA516uL)),
                ((-1, 8, 0xAEF226F534FE1CBBuL, 0xEF4411E16E67BEA3uL), (-1, 7, 0x875791395537ECE2uL, 0x7E1A4E8AE013DE1FuL)),
                ((+1, 8, 0xC33B5043CBE547D8uL, 0x583F054B106ADF43uL), (-1, 8, 0xC70A272B38DD8DB0uL, 0x524B16D52A7701C5uL)),
                ((+1, 9, 0xFC13DCCB47C8F08AuL, 0xD017A1C03F1B53C6uL), (+1, -1, 0xF9EAFA364522E58CuL, 0x68EBF03DF4560A48uL)),
                ((-1, 8, 0xA5BE0A8383D79891uL, 0x5DE3D797C8532ACAuL), (+1, 9, 0x831E16024AD76680uL, 0x94880510BA94251BuL)),
                ((-1, 9, 0xE91E373747AF2B7AuL, 0xA06FFECA4C1A1835uL), (+1, 6, 0xD7AD00BAE8297A2EuL, 0x2715639AEB4CFCE3uL)),
                ((+1, 7, 0xA0ECA4B707C537A5uL, 0xDA9611012CA1E3D4uL), (-1, 7, 0xE5A24152C73923BFuL, 0x05D342BFF593B038uL)),
                ((+1, 7, 0xD71916CF8F37B64EuL, 0x8534CD1B406AC9B2uL), (-1, 4, 0xD478018A79E4C07CuL, 0x422B233959E945CBuL)),
                ((-1, 4, 0xAA8791987ED14A7EuL, 0x6FF42546CAF78342uL), (+1, 4, 0x80DC416A9A8D4322uL, 0x1F8DAC4C99F4B483uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p375_0p4375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xA82069F21FE35E50uL, 0x212BC28C19205FF2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0x8C875BDA6AD0BB33uL, 0x7AF1F388D5624430uL), (+1, 2, 0x95F501B079CB8417uL, 0x1821EFC7DEE2EA77uL)),
                ((-1, 4, 0x9BE45E4774CC255EuL, 0x1AFEC13DFA34515DuL), (+1, -3, 0x82C57D5FCD8C0D21uL, 0x9D64D910487089D2uL)),
                ((-1, -1, 0xCEBB747AE814DFD6uL, 0xA818966C58F87BDEuL), (-1, 4, 0xB663727F9CA33201uL, 0xBF587245F2E60118uL)),
                ((+1, 6, 0x8A4447F6943A61B4uL, 0xB2C4F24C33356F1FuL), (-1, 3, 0xDE1BD8526F9EAD11uL, 0x877D649BA2595264uL)),
                ((+1, 3, 0x984F39712CF6D119uL, 0x15F7DDD0B4D8D88AuL), (+1, 5, 0x8F056B61EFFF57BCuL, 0xD299AB5A9A5938D6uL)),
                ((-1, 6, 0xA50C04725C485BC3uL, 0xCF461ADA4D80DD8EuL), (+1, 4, 0xAC96B31C887F2732uL, 0x3DCCB452FE7F571FuL)),
                ((-1, 2, 0xC0064D187C3694A8uL, 0x5A17D0EFAA2A1172uL), (-1, 4, 0x96B0B703430248A1uL, 0x16C18737C321503FuL)),
                ((+1, 4, 0xCD06BB8C3EACF7B2uL, 0x97C90053D2A8FAE5uL), (-1, 2, 0xD775916D4C046276uL, 0x047C0164DF0DA952uL)),
                (Zero, (+1, 0, 0xC9C64DB57220B677uL, 0xC6CC6FB4AF267F73uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p4375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xE888E31EC9D30D6BuL, 0x63F60346A6781E5AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 1, 0xD16FF17B14818591uL, 0x55668683038267D0uL), (+1, 1, 0xB37F3AF255828603uL, 0x9FD6AED6BFD5A9CFuL)),
                ((-1, 3, 0xCEE84E6D09D8A463uL, 0xE637CA6CAC0CC8F0uL), (-1, 2, 0xC9096EF069BA5459uL, 0xDDC73A0B4C8A761CuL)),
                ((+1, 4, 0xBD12DDCE541468BFuL, 0x1501DCF478916A27uL), (-1, 4, 0x9B50E83EA495B3D4uL, 0x0A829929E01C57EBuL)),
                ((+1, 5, 0xD47184699426D0EEuL, 0x6B98A1BE9FF79116uL), (+1, 3, 0xD019818DF4FAD744uL, 0x9849B0A20C0C5CD3uL)),
                ((-1, 5, 0xE8BED20DD6D61D34uL, 0x103AD5023308CEB5uL), (+1, 5, 0xB2D8DF9C734DA5C7uL, 0x1182E588345D948BuL)),
                ((-1, 6, 0xA252D5FF67E16577uL, 0x5482949F8C26E24DuL), (-1, 3, 0xA10777C4BB80784DuL, 0x21F608C385B70451uL)),
                ((+1, 5, 0xDA18BDC73AAC20BBuL, 0xF23D2072B73A2BE1uL), (-1, 5, 0x97C8C430670D649AuL, 0x77C76D2D86B7DC8EuL)),
                ((+1, 5, 0xA3E97D1ECA477A71uL, 0xE02E9318DBA09A5DuL), (+1, 1, 0x94E4E009722FA145uL, 0x9FF9999846C87EB5uL)),
                ((-1, 3, 0xE3922F97584F6264uL, 0xA0B0417DD3039C5BuL), (+1, 3, 0x89792BA37CDEE661uL, 0x546A16D2B200F42CuL)),
                ((-1, 1, 0xCFE1CF58906BD182uL, 0xC0A0FEB7E965127AuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xB5CFBD01F2E52D51uL, 0x63E7A633C7038B6FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xE248829D9987CFFBuL, 0x19E6AE2E6CAE713AuL), (+1, 0, 0x9C225E8C3E8E4BAFuL, 0xDF157F05A7083507uL)),
                ((+1, -1, 0x81506C063A2F29CEuL, 0x8E2677B5ADEE478BuL), (+1, -1, 0xBA6ED8FA50F156B4uL, 0xA5DD10D1E2CEF85BuL)),
                ((+1, -3, 0xC83C07EE8AA636EEuL, 0xA5AC4C3A8FE52BA7uL), (+1, -2, 0x96060D9F03DEED61uL, 0xA678D3B84E086183uL)),
                ((+1, -5, 0xE8EC773ED906A65EuL, 0xD736DAAC6B7273A3uL), (+1, -4, 0xB1DA005737B4A83FuL, 0x9F8B0A93444DB4EFuL)),
                ((+1, -7, 0xCA9D396E7CCAD3DDuL, 0x53F7B09AA163EFE0uL), (+1, -6, 0x9F786E5A4364A208uL, 0x1B45DF842D9FDBC4uL)),
                ((+1, -9, 0x8B3A014F181D76D3uL, 0x6B55A1D04DA5D8BDuL), (+1, -9, 0xDC87C8C67E61EC0FuL, 0x8BB72D0C418694A7uL)),
                ((+1, -12, 0x9418D876DF60C703uL, 0x148F57D61EBC6125uL), (+1, -12, 0xE9086320F5228C41uL, 0x67F482D3E6B389BBuL)),
                ((+1, -16, 0xE49A4C9FC09B0EBFuL, 0xFA41BDCD82771A86uL), (+1, -15, 0xB756E216BF58845AuL, 0xCFFCA26A804ABCEEuL)),
                ((+1, -19, 0x918D9E3E2CE912FBuL, 0xD43BBBADAD3C7A7AuL), (+1, -19, 0xD09A6011FC355153uL, 0x3CBAC9CCA4697A83uL)),
                ((+1, -24, 0x9ED872222F9BEA2FuL, 0x79CAECC570A2630AuL), (+1, -23, 0xA0674684F2DA6D3BuL, 0xADB1C1D4940CD69BuL)),
                ((+1, -28, 0x83B06B07F89EB895uL, 0x597E214C3A624D2CuL), (+1, -28, 0x92E0D031F6717B8DuL, 0x05DC114100E8F51EuL)),
                ((-1, -36, 0xD4B89CD0AA8DF8F7uL, 0x095C0D651355D84DuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xB4C6136EFFF0E299uL, 0xD3CD55A5EC532678uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xE9F98665F876B379uL, 0xB42C33B1BFE3DB68uL), (+1, 0, 0xA8B6C16F45472E16uL, 0x9F85E7281E35E0ACuL)),
                ((+1, -1, 0xA02A3A61A83876E2uL, 0xFE684EE7E1CE77C4uL), (+1, -1, 0xED154608E39229C8uL, 0xBF38CA8606A5A15FuL)),
                ((+1, -2, 0x95E16988F7C2066DuL, 0x8214B78B095FD0A2uL), (+1, -2, 0xE210EDFE68FC1706uL, 0x8BBD8D3C0989E522uL)),
                ((+1, -4, 0xD00B63F761E9F135uL, 0x21CEC3852869BF71uL), (+1, -3, 0x9FA48EAA7117CF6CuL, 0xD8F43C5B0F803551uL)),
                ((+1, -6, 0xE0CCCD893076ECC8uL, 0x3E7909B525ADB093uL), (+1, -5, 0xAEB5538C9EE71933uL, 0xC8314FD619528073uL)),
                ((+1, -8, 0xC1468A45AF414D14uL, 0x59F9982E4FF8BCE9uL), (+1, -7, 0x9781405E08741EABuL, 0xFDB472BD9E40022CuL)),
                ((+1, -10, 0x85D6D527A750B2EDuL, 0xBCFB5FA091C22150uL), (+1, -10, 0xD2C766995081245FuL, 0xF75FB021EA48B484uL)),
                ((+1, -13, 0x965A8975E202B6FEuL, 0x9B9D3E1C63C574B6uL), (+1, -13, 0xEC903FDD9435EF51uL, 0xAEC2374BC263AEB9uL)),
                ((+1, -16, 0x888CA7EE7178002FuL, 0x4C846C77F740508CuL), (+1, -16, 0xD695C5BBC6546922uL, 0xE94CE561D952F3C4uL)),
                ((+1, -20, 0xC95091733107CBF3uL, 0xB3CA846D7706CA7AuL), (+1, -19, 0x9D8A799ED92A65E1uL, 0x2C20ED1D268F2179uL)),
                ((+1, -24, 0xEDC7BB7B4F3C1D2BuL, 0x3FFDDF3EBA8C94B7uL), (+1, -23, 0xBB3B8D4C9F82768CuL, 0xFAD3A00D5F5507DFuL)),
                ((+1, -28, 0xE42D426349FEEBB9uL, 0x8832EB95AD06AB7CuL), (+1, -27, 0xB395CC804F45D10CuL, 0x827FA84D03F957A8uL)),
                ((+1, -32, 0xAFCA05A2E35BB324uL, 0x81043D9873B22558uL), (+1, -31, 0x89230FEF0266345FuL, 0x205814031FC11F37uL)),
                ((+1, -37, 0xCF27CB41D12DA314uL, 0x203D59F0F47A2E06uL), (+1, -36, 0xA41D57C8C5594F9FuL, 0x1FFC55224B7F3EE0uL)),
                ((+1, -42, 0xB44511D3A12ACC2AuL, 0xCF2DCB098E83F53FuL), (+1, -41, 0x8C891B85D5935300uL, 0x92AA6FCEFD91CE3BuL)),
                ((+1, -48, 0xD3A383773067E6B4uL, 0x56C6ACD13990E51BuL), (+1, -47, 0xA6DDAD6117764408uL, 0xF4A282CF03C8577DuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm8_12 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA5F11067D6600A2EuL, 0x1049BC06BEF18C3CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x86A05A318B901010uL, 0x66284828B100A8C7uL), (+1, -1, 0xD22CB019851BFFF5uL, 0x026F931073482C52uL)),
                ((+1, -3, 0xCE1E613A758E441DuL, 0xA7581EF731DA4A49uL), (+1, -2, 0xA1D47DD3F457C7A0uL, 0x24ACB794E457922BuL)),
                ((+1, -5, 0xC54413A9298806BAuL, 0x17F5AD4A3E986995uL), (+1, -4, 0x9B1041164BD9F508uL, 0x418D29A93A2F54F4uL)),
                ((+1, -7, 0x83E54154C5CE1583uL, 0xCA470D1128557AC3uL), (+1, -7, 0xCF360B2833456AA3uL, 0xAB6B803B9FA8C121uL)),
                ((+1, -10, 0x82BDAC3B4C00CD89uL, 0xE015AB5CDAA5178EuL), (+1, -10, 0xCD4C22CF44393B8BuL, 0x5B676D468A0C7187uL)),
                ((+1, -14, 0xC7437D2C69BDC880uL, 0x161095A1681A3A1BuL), (+1, -13, 0x9C80560427FD0527uL, 0xA28A1870F84013AFuL)),
                ((+1, -18, 0xEF0279292BB9A424uL, 0xCE3518EBA7A345B1uL), (+1, -17, 0xBBC46F2B3F44A243uL, 0x08287C33F838057DuL)),
                ((+1, -22, 0xE46BC173683B2C89uL, 0x4675F2EED2BEC85EuL), (+1, -21, 0xB35D492D807BB63DuL, 0xB37D27FA5FE19845uL)),
                ((+1, -26, 0xADF255BEA0886DE5uL, 0xE36B08A885A4CB67uL), (+1, -25, 0x889ACEDA9F86AFB7uL, 0x7818A8E78C69016BuL)),
                ((+1, -31, 0xD0C48F6C9EBDC974uL, 0x652EE0AA03104B6BuL), (+1, -30, 0xA409398A7B6802DEuL, 0x0B6AD6E0A4513D2DuL)),
                ((+1, -36, 0xBFCC160E33064397uL, 0x14DDC34D4F9C8D6FuL), (+1, -35, 0x968AE4E369138BFCuL, 0xDD8F115FB976DCF2uL)),
                ((+1, -42, 0xF911850B1AF03BB2uL, 0xE04D557714C29FAEuL), (+1, -41, 0xC3BDA05FE98FB881uL, 0xBF067AD8C1257452uL)),
                ((+1, -48, 0xBE445ADF3E10CDE3uL, 0xAA8098C0A9764912uL), (+1, -47, 0x955E5E2A5357675CuL, 0x8312CF6DBBDF140FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm12_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA345531D71DD96DAuL, 0xBAC0036A7A51F5C9uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x9D6D88A16B62609DuL, 0x8961A2934A4E1BCCuL), (+1, -2, 0xF7654F601ED5073EuL, 0x9D6F37A35519125FuL)),
                ((+1, -4, 0x8F886D0A9924540BuL, 0x07D3BC5E81C26864uL), (+1, -4, 0xE17A00047EEDCF8CuL, 0xA24B35456FED034CuL)),
                ((+1, -7, 0xA49C758F99B96232uL, 0xF918DA77E998FFE5uL), (+1, -6, 0x814363B003EA9948uL, 0xC9473074401B7643uL)),
                ((+1, -10, 0x8517D0AB1E331FEFuL, 0x8E22A5D49BEEC2E4uL), (+1, -10, 0xD1159B272C74AAAAuL, 0xB02D886F989D593AuL)),
                ((+1, -14, 0xA0DE4D6DE3A9FCFCuL, 0xE3524E191101AAA0uL), (+1, -14, 0xFCB0DD14F7CE009FuL, 0x087E7C32859534F9uL)),
                ((+1, -18, 0x95AC6E4749632F79uL, 0x8CAFBC7DFB21B3E4uL), (+1, -18, 0xEB163B44D9CEE316uL, 0x58A7BF7F2755E0F7uL)),
                ((+1, -23, 0xD8C00B802D4846B9uL, 0x6B2B0948BFE0E327uL), (+1, -22, 0xAA41F278B276CDD4uL, 0x79B7FA96FF3BE961uL)),
                ((+1, -28, 0xF2EEEA35254D62DCuL, 0xE444BF5DE7EF53A4uL), (+1, -27, 0xBEC5CF7EB74A3DF2uL, 0x3B443267FEA48D7BuL)),
                ((+1, -33, 0xCC828C8FABDF6E0DuL, 0x1D2CC66AD6FA7307uL), (+1, -32, 0xA0A40E876A2B8557uL, 0xD456158411044DBEuL)),
                ((+1, -39, 0xEF87C9844BF328E7uL, 0xAA7F8ACF6D301594uL), (+1, -38, 0xBC1C67A0192BF437uL, 0x7B3F4FFAC97DAB9CuL)),
                ((+1, -45, 0x9F65B3D0BD611798uL, 0x7712666A8F2049E0uL), (+1, -45, 0xFA6471AEF281074CuL, 0x373921E07F5D6CD3uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA300045F188C308BuL, 0xC9C9BB8B8C8C0A79uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x9A62E786156A7470uL, 0xED3A72165B6ED957uL), (+1, -2, 0xF28587EA7B1DB2B2uL, 0x254BBD66216C2EEEuL)),
                ((+1, -4, 0x8AF71C99775A3EA3uL, 0xA1ED986ECDF483C1uL), (+1, -4, 0xDA48F614E8DDE1CEuL, 0xEE643BFC2194D480uL)),
                ((+1, -7, 0x9F05845505891B97uL, 0xF5A09B310AFF2F35uL), (+1, -7, 0xF9C9C443D9ED7EDCuL, 0x34F8E78EA4B37766uL)),
                ((+1, -10, 0x824368B8D294C06FuL, 0xEBEB2EBF16B47E6EuL), (+1, -10, 0xCC9E80D2A0EAD198uL, 0x3BFB620856EECB60uL)),
                ((+1, -14, 0xA2CC55DCB847FCECuL, 0x766E227BCB2217CBuL), (+1, -14, 0xFFB8B5A194D68032uL, 0x34FECAF98150408DuL)),
                ((+1, -18, 0xA139EF9EC4ED458DuL, 0x6829057F8A1A4897uL), (+1, -18, 0xFD40E23C6E4973A9uL, 0x01461E83874C8554uL)),
                ((+1, -22, 0x8190B6ABC7FCAE1DuL, 0xDD979E4539C0A303uL), (+1, -22, 0xCB859A287D4811CAuL, 0xCE89F26FC4846C2EuL)),
                ((+1, -27, 0xAB6453D552EFBDB1uL, 0xDBA13674552D8175uL), (+1, -26, 0x869C12847B321F10uL, 0x353C7CFAF25D9818uL)),
                ((+1, -32, 0xBBD5CCAB3A97BF45uL, 0x47E02C97EF2C6390uL), (+1, -31, 0x9386E5C625B9C8F6uL, 0xE7109BEB76CD289CuL)),
                ((+1, -37, 0xAA7D85AB27C2C96CuL, 0xC462211332FDCA70uL), (+1, -36, 0x85E6C05523E8734EuL, 0x3C69DD70A6D24493uL)),
                ((+1, -43, 0xFE2250004B75456AuL, 0x8A9FFCC8B494E19BuL), (+1, -42, 0xC7992ADA2A021D8FuL, 0x36E87C24AA267F63uL)),
                ((+1, -48, 0x986FE9775245003AuL, 0x921912EB6EEBDBB1uL), (+1, -48, 0xEF721FC9A60A4A5BuL, 0xFD3A63B974F25FF1uL)),
                ((+1, -54, 0x8D700C4E05719200uL, 0x1DC071E42541B53BuL), (+1, -54, 0xDE2BD5BFB787018CuL, 0xFF88F67756ECAF85uL)),
                ((+1, -61, 0xBA51869C70779DD7uL, 0x5E8302E1E6787FF2uL), (+1, -60, 0x92555C0613199249uL, 0x608E30566D9F99B3uL)),
                ((+1, -68, 0x8C7CBDE51AEEE981uL, 0x0216FFAFC33E01E2uL), (+1, -68, 0xDCAD880006F1156FuL, 0x3E2A5EB84A2AE76EuL)),
                ((+1, -95, 0xC7291DC616556C09uL, 0x306685A1EB61BC23uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9837BDEA8F25BuL, 0xECA814B3F4F458E5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF47FC81BCE4EF564uL, 0xCBF630D9217243AAuL), (+1, -2, 0xC007787662E50A90uL, 0x7BE3C4C884228449uL)),
                ((+1, -5, 0xB4372186553EB724uL, 0x3B357E5BE1679547uL), (+1, -4, 0x8D8A726FDE17640AuL, 0x9ACDD2FDA4D2F164uL)),
                ((+1, -8, 0xAD7E86D2B0286F9CuL, 0x34502857AAF87C3BuL), (+1, -7, 0x884316759098C7D7uL, 0x427AA10E95A0D225uL)),
                ((+1, -12, 0xF49DF57C3B44F888uL, 0x3764ED41592AE655uL), (+1, -11, 0xC01F2BE4CBEF41C8uL, 0x0954FF1FAF90FE13uL)),
                ((+1, -15, 0x8635FA454619A189uL, 0x77522AEB5BC638F3uL), (+1, -15, 0xD2D162956F807F00uL, 0xB1946464652F8E34uL)),
                ((+1, -20, 0xEDCA5926CB983FAFuL, 0x1D54119DA7095594uL), (+1, -19, 0xBAC299ED9EDF6C3BuL, 0x6DC007401EE69F86uL)),
                ((+1, -24, 0xAE101EE36BF23991uL, 0x6FF6B7726328C327uL), (+1, -23, 0x88B56FDB3FDC6053uL, 0x4AF5D03FCC9DAD34uL)),
                ((+1, -29, 0xD5A98AB9956A89ABuL, 0x82B72518258F9722uL), (+1, -28, 0xA7CF5953DADA62E8uL, 0xA78B57845235A66AuL)),
                ((+1, -34, 0xDDB76666B1948E16uL, 0xCE42EC9E74AB10C4uL), (+1, -33, 0xAE22BA80A1C9F4ECuL, 0xDC94B58C028C9FA8uL)),
                ((+1, -39, 0xC305704245F66206uL, 0x890E5606A1A3E59EuL), (+1, -38, 0x992B59038F13DA2DuL, 0xA5E8702C4AB9796EuL)),
                ((+1, -44, 0x90F713C73687DF7EuL, 0xA48F5D7FABE27B90uL), (+1, -44, 0xE3B5F19FB7BCA46CuL, 0x7A097F20C99FF185uL)),
                ((+1, -50, 0xB44A3EB6A6389F45uL, 0x5105C7F10BFA8226uL), (+1, -49, 0x8D997590719ABEA7uL, 0xF3735EADC6480F92uL)),
                ((+1, -56, 0xB7B9FB08C88BFF11uL, 0x543C06C488A2FF87uL), (+1, -55, 0x904C66E85B9BB5D9uL, 0x39E5231428F8DC20uL)),
                ((+1, -62, 0x9371283CFE323B86uL, 0xE59EB612340D9C27uL), (+1, -62, 0xE799F449F537F0F4uL, 0x870531A47DF4AB49uL)),
                ((+1, -69, 0xAB86E546EED66734uL, 0x55A3F79141DAEE68uL), (+1, -68, 0x86B78969A2595E2BuL, 0x87F8ADA88D072304uL)),
                ((+1, -77, 0xE9259DC32B48CF83uL, 0x85016609BDB871C5uL), (+1, -76, 0xB71CF9302BBBD597uL, 0x6560540D0479E803uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E441545uL, 0xAB70FCC4E1436397uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8EA9AAAE1BC02062uL, 0x1472A9B1639BBFEEuL), (+1, -2, 0xE018199512A3FF60uL, 0x828EC551CC6748F3uL)),
                ((+1, -5, 0xF5EE2BBCCF7E9438uL, 0xCEC064AD66FA99ADuL), (+1, -4, 0xC1273B54610957E5uL, 0x3F753B1340BD0A80uL)),
                ((+1, -7, 0x8AEC4E54A39C792EuL, 0x697775058A4660DBuL), (+1, -7, 0xDA3847FC8394CB64uL, 0xDF3E3DF3D0548DF2uL)),
                ((+1, -11, 0xE702003CA993986BuL, 0xE85DAE74E534801FuL), (+1, -10, 0xB56EE097A6506B32uL, 0x95E7A4B38FD0FCD6uL)),
                ((+1, -14, 0x966F69F4DF64233EuL, 0x1C5D930A2EA0B8D5uL), (+1, -14, 0xEC4D96653E386B8CuL, 0x487D3A90E4C6AEF5uL)),
                ((+1, -18, 0x9F6F57821F10457AuL, 0xFB2AC7E2ECF81CBDuL), (+1, -18, 0xFA7096CA2998BF86uL, 0xDB675DB4310BB4BFuL)),
                ((+1, -22, 0x8CFA5B42AE74AF87uL, 0xE2609A43C6937F8AuL), (+1, -22, 0xDD7299725CD8ABAAuL, 0x940268723C3EED89uL)),
                ((+1, -27, 0xD3747404DC0A83D5uL, 0x0A4DB212A2CBA760uL), (+1, -26, 0xA6138781B3233F41uL, 0xBC304195CD825BCBuL)),
                ((+1, -31, 0x881A739C28A6C8D8uL, 0x65883F20FB0605AFuL), (+1, -31, 0xD5CA651E3609279BuL, 0x446E5443E9BADB63uL)),
                ((+1, -36, 0x965FB7B90336AB2BuL, 0xCDE69CD82F504357uL), (+1, -36, 0xEC34EE918C05981EuL, 0x0CDADE34AC204E53uL)),
                ((+1, -41, 0x9245C107166D3882uL, 0x704D6E116EF551B3uL), (+1, -41, 0xE5C3A721D3423512uL, 0x89EB579E95BDB77FuL)),
                ((+1, -47, 0xE3AC1CD1697F934CuL, 0xCB21212D6D0FEE7EuL), (+1, -46, 0xB2D03C26DD4543D6uL, 0x62881214C74F4FBCuL)),
                ((+1, -52, 0xBA0A9690322B68F2uL, 0x0E27BB3F79DE1469uL), (+1, -51, 0x921DD5B8E10896C2uL, 0x0A71E529CC2A7D80uL)),
                ((+1, -58, 0x869BD76D945F182CuL, 0x311D6A4407B0AA0CuL), (+1, -58, 0xD37164793D1A56E0uL, 0xF1DCAD2E13279FC9uL)),
                ((+1, -63, 0x88F94EC970182D9FuL, 0x9BFEC5493189D087uL), (+1, -63, 0xD72874E5A7184507uL, 0x4DCCA298D339789FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -1, 0xE08A3FD6E0EE88CCuL, 0xB0AE856BC5DE6DE9uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 4, 0xE009C65D7D172A30uL, 0x05DDFE193A3EFD72uL), (+1, 5, 0x9235E2C82A009E67uL, 0x34241D2D25760FB6uL)),
                ((-1, 8, 0xAAF41378F43F03F6uL, 0x5FAC2BF0BABEE5E8uL), (+1, 9, 0x8A08942FBD5409C9uL, 0xA069642BEBB77ED5uL)),
                ((-1, 10, 0xE095C0B586F8F525uL, 0x7512CD682F317D15uL), (+1, 12, 0x897080D264315C00uL, 0x4A1AE764CDE37BECuL)),
                ((-1, 10, 0x8D4D726D86289793uL, 0xA7B8F073FF8FDE51uL), (+1, 14, 0x962C164F39B6EC61uL, 0x554B0D8ED93789D1uL)),
                ((+1, 14, 0xF020DF47613FCDD7uL, 0x0B7F2F6364B0D18EuL), (+1, 15, 0x9E85DD081B3E8A6DuL, 0x7F96A40F9E6764B7uL)),
                ((+1, 17, 0x82D26B1C7CF6FC36uL, 0x2BF120C246B4E6D6uL), (+1, 12, 0xF36A068C168AF3BBuL, 0x6741E31A50CA422CuL)),
                ((+1, 17, 0xA9B78BC664A6A8B3uL, 0xD5B6B6CCA1EB82E2uL), (-1, 16, 0xE54CFB6EA032E1E1uL, 0xE3A07A65A108F63DuL)),
                ((-1, 16, 0xF011BCC3BC0788ACuL, 0x542D21DBD4FB1ADAuL), (-1, 16, 0xF7000BE4346BCBCCuL, 0xA253DA8FC2070FAAuL)),
                ((-1, 18, 0xBA0A64E6F17A88F3uL, 0x0B9D62F3EC2ACEA9uL), (+1, 16, 0xC2FA7F49FE199B8DuL, 0x5DDE6A33DFF9254FuL)),
                ((-1, 14, 0xAD6978FCB243E44AuL, 0xDDCED6C6C3AE4799uL), (+1, 17, 0x828F0A3548D9A15FuL, 0xCFFD1AE64065972CuL)),
                ((+1, 17, 0xCAC44B2CD8E0526BuL, 0xA4CDFE49382563A5uL), (-1, 15, 0x9FE79B9E8EEF63BEuL, 0x0D16F79AFA1E2D6DuL)),
                ((+1, 13, 0x928E232D01D01E60uL, 0x4F442CA90C394711uL), (-1, 14, 0xE66F07478C006A7CuL, 0x69D885A832E58A59uL)),
                ((-1, 14, 0x8B0CFA1E862A14D8uL, 0x96CB833E58ACCBDBuL), (+1, 12, 0x8D866EB90A8D7C10uL, 0x8CC50DD312F48AF7uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -2, 0xD5E5435E09BA5BCEuL, 0xFD7DF3546946D795uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xAE9FE92F7F395D13uL, 0x71A13CE405A6A6BAuL), (+1, 3, 0xBB3914E19B44437EuL, 0x97D1461D1DB3BF4AuL)),
                ((+1, 4, 0xB87A544D75DDC4CCuL, 0xA45920D22A40B073uL), (+1, 5, 0xAEC9608954FF3763uL, 0x0A328CAD152E3CD6uL)),
                ((+1, 7, 0x97F653667767F11AuL, 0x24F68C999C3E80EEuL), (+1, 4, 0xCBC419387CA816CDuL, 0x8FBB15B880631F77uL)),
                ((+1, 7, 0xF490F33CFC5CC1DCuL, 0x41DEAC9428A7C6A4uL), (-1, 7, 0xAC80D8614CAB0D98uL, 0xD88870C37A7F29C9uL)),
                ((-1, 7, 0xF86DB3E4E68FD69AuL, 0x5E42BB37461E8407uL), (-1, 7, 0xF0F5E84133EB6842uL, 0xB60A7F17245DBEAEuL)),
                ((-1, 9, 0xCC8D0DB0DC009807uL, 0xD7D95E8FD11E64ADuL), (+1, 7, 0xE5A1B4A7ADDD8A08uL, 0x32CC04CD4B9B51B4uL)),
                ((+1, -1, 0x8A8FA73A30F8A0C7uL, 0x8CC85CFAFE156969uL), (+1, 8, 0xBFABB1CE07397CF4uL, 0xF5C76B8E0EF5AC80uL)),
                ((+1, 9, 0xBDAD7922ACD34EFBuL, 0xF854CED9B096DABDuL), (-1, 7, 0x945F032BF8AF0B4FuL, 0xBE5FFDC3D0800A32uL)),
                ((+1, 5, 0xF04830F62A61916AuL, 0xDCD8CF0E87816456uL), (-1, 7, 0xB004141E486C39D9uL, 0x63D7C656B10E6764uL)),
                ((-1, 7, 0xB8F4160F344A1D16uL, 0x75C91D8E570F59CDuL), (+1, 5, 0xA262C1980088258EuL, 0xD24725EDE8CAA64DuL)),
                (Zero, (+1, 3, 0x9F9402B4E79A6642uL, 0xD37650089F16C626uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0x996B31177B2D3822uL, 0xDB6BF648ED48C8B2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0xFD7CA6F8FB7F04C4uL, 0x98CD1CDF7925BB86uL), (+1, 1, 0xA93F33FAB1487993uL, 0xA233E18A5BF98383uL)),
                ((+1, 3, 0xC744824C16537E05uL, 0x6CF5AF40963EA8AFuL), (-1, 3, 0x8F72A5D7CEDF9A1CuL, 0x869EE69DFE56849EuL)),
                ((-1, 4, 0xADDCCA1E11089B87uL, 0x3ED790037F4364CCuL), (-1, 4, 0xA8539057851BA337uL, 0x8A945CCEF4F42092uL)),
                ((-1, 6, 0xA0411370387EE66DuL, 0xD2EDFCB3CCC35E58uL), (+1, 5, 0x88A48F442CDDA517uL, 0xAEFED27127DF931BuL)),
                ((+1, 5, 0xC5825225176D3F5FuL, 0x7EBB1D60FD6E022FuL), (+1, 5, 0xE55B4CDCCE7F3119uL, 0xA6979C2472C19D55uL)),
                ((+1, 7, 0xADE2A31F8E227A01uL, 0xF07566976CF3F167uL), (-1, 6, 0x8446094D694CA875uL, 0x9D3538D07A152994uL)),
                ((-1, 6, 0x852C7EEBADCD8502uL, 0x84B0F8773ED955C7uL), (-1, 5, 0xE9F6CC0BBF60754DuL, 0xCF8124C80CB9861AuL)),
                ((-1, 7, 0x89AE1DE83848947AuL, 0xBB77BB4288948545uL), (+1, 5, 0xE19BC58F4D1DBB8CuL, 0xB263A6513C024788uL)),
                ((+1, 5, 0xAA8950114AA9ACBAuL, 0x03C36F96EDED31EAuL), (+1, 3, 0xFF59585DE6CF40AFuL, 0x53B8142B8D9768A5uL)),
                ((+1, 4, 0xD7264133EC85485AuL, 0x112972BFD563CF44uL), (-1, 3, 0xD157567B5CBC6C24uL, 0xC9F02810EFB43FC0uL)),
                ((-1, 2, 0x971A2545AF7B6F54uL, 0xA0A5547FD8C6E1E7uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -1, 0xE08A3FD6E0EE88CCuL, 0xB0AE856BC7209693uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xD159AC31CBF42D15uL, 0x6F8AF1CDCA2643E9uL), (+1, 0, 0xBB660D1C0CEC7C5AuL, 0xF1871DC76DA48863uL)),
                ((-1, 0, 0x9C1C65045D2C67F6uL, 0x1F5D5F2372F53698uL), (+1, -1, 0xE478C15021E762A1uL, 0xDAD71537DA797993uL)),
                ((-1, -2, 0xF4FBDE9285B14B4BuL, 0xA9F332E60541406CuL), (+1, -2, 0x963A5005B8F05E20uL, 0x05C6E627689EFA98uL)),
                ((-1, -4, 0xDD286960DB7DD3E4uL, 0xF56AC066D497F4EDuL), (+1, -5, 0xE6ADCCDA45437BABuL, 0x630392208A89FE74uL)),
                ((-1, -7, 0xEA13346123FC41CBuL, 0xDE56EFF6A2279F09uL), (+1, -8, 0xD11A43F2B962819CuL, 0xB278C1025667EC6DuL)),
                ((-1, -10, 0x8CEE4CA331F3DF05uL, 0xCFF6B4B2CCE397F6uL), (+1, -12, 0xD77C7227F566F286uL, 0x56E1899F802B7A3FuL)),
                ((-1, -15, 0xB0EACCF6545C2393uL, 0xAE28DA9B02B04F85uL), (+1, -17, 0xE4EFBFD4F4B05B23uL, 0xEA4410875D52ECAEuL)),
                ((-1, -21, 0xBF2DAA4986855484uL, 0xEA21665AF160B665uL), (+1, -23, 0xCAE231CA56344001uL, 0x47961F18642543ECuL)),
                ((-1, -29, 0xD9FDBE03C35E16B6uL, 0x23D7EC9FC8F5FE19uL), (+1, -31, 0xABE899E95454AAC8uL, 0x1536824C6B1A4EC5uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 0, 0x9569334F8B508630uL, 0x2B6DF06F000FD27CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xEC3F44E730A8B6BCuL, 0x147A99FBFB41AF0DuL), (+1, 0, 0xAFE4E881F55B763BuL, 0xCCCDBBE0F0FDB6DDuL)),
                ((-1, 0, 0xA01853DF446A999CuL, 0xAD4FB833F273F5E6uL), (+1, -1, 0xD1F2BB1A1018D3F1uL, 0x31841826D187E4FFuL)),
                ((-1, -2, 0xF4EFE4B72A54970EuL, 0xB3163AFD81F00C2BuL), (+1, -2, 0x8ECBC6E36ADF53F3uL, 0xD475A016056E98F9uL)),
                ((-1, -4, 0xE9FE387D7C9F4597uL, 0x0A3DA8BED026DD8FuL), (+1, -5, 0xF43100FB319C3237uL, 0xA09B90E7B230C591uL)),
                ((-1, -6, 0x91A1DD5EE90F2F39uL, 0x834C68B3314CC354uL), (+1, -7, 0x8898427FA6EC92EFuL, 0x1550AAD0C5376D55uL)),
                ((-1, -10, 0xEEFEA6065FD4FDF5uL, 0x3744FD1583D979A9uL), (+1, -11, 0xC9D93238968536D9uL, 0x94F4F4677360F1DEuL)),
                ((-1, -14, 0xFFCB5D21904820BEuL, 0x37B13D9DDD1452FBuL), (+1, -15, 0xC26AB85916F870A9uL, 0x85C975FEED90C7C5uL)),
                ((-1, -18, 0xACBF2FC59A8E199DuL, 0x4245E11C4CF0FFFAuL), (+1, -20, 0xEB765540C77AF5ABuL, 0xB2175B3F95BA5092uL)),
                ((-1, -23, 0x8A5799358E135481uL, 0x97BCDDCFA7AB17ADuL), (+1, -25, 0xA7BD602AAAC17041uL, 0x72C9058568B2845DuL)),
                ((-1, -30, 0xEB5B4F1308965105uL, 0x030A9D23E10A4776uL), (+1, -32, 0xF9EA3E8589791D5FuL, 0x156D3A1CA11F153FuL)),
                ((-1, -37, 0xAC5F800AC2FBDAA3uL, 0x5B66A0BE49C374EAuL), (+1, -39, 0x9B21466EEB1B2911uL, 0xEF8F667D0A8DC31EuL)),
                ((-1, -46, 0x82AFBCEC8E2A07FAuL, 0xCF817521FA55BD38uL), (+1, -49, 0xB5A5FC0D8D61F0D1uL, 0xE2ADE3DD2F859014uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 0, 0xE44915D1249ADF91uL, 0x684CB40E700C6DEAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xAA14258BCA289172uL, 0x868DB5EE6D51F670uL), (+1, -1, 0xAFED2AE845139D06uL, 0x39A712C5DBEAFBFEuL)),
                ((-1, -2, 0xDBA21160E82F6282uL, 0x9F5DDF9ABE9F2C23uL), (+1, -3, 0xD1FE72439F2F1993uL, 0x85B92037ED392D36uL)),
                ((-1, -4, 0xA1665A49132D8DFCuL, 0xC7986A531BEC1C21uL), (+1, -5, 0x8ECED3A8CEC2AD58uL, 0x9F11C75DC92A8945uL)),
                ((-1, -7, 0x950A0091188ED8E0uL, 0xCFF6D74347FBBDD5uL), (+1, -9, 0xF421FC51747635B2uL, 0x1C885DCB5074A68FuL)),
                ((-1, -11, 0xB42D1A141EAD27B5uL, 0x00427CB2AA38903AuL), (+1, -12, 0x887D3EAFF3881840uL, 0x31509A2DB26E6BB7uL)),
                ((-1, -15, 0x901CABF2888386C9uL, 0x6D4D39954959D001uL), (+1, -17, 0xC98A081EDE8C9E74uL, 0x23C62BEA18316AC6uL)),
                ((-1, -20, 0x96C55B42A06722E0uL, 0x8E2947D33786A0D6uL), (+1, -22, 0xC1EC8685A87F7E90uL, 0x136F9103A3E27A79uL)),
                ((-1, -26, 0xC777D658E61E3FEAuL, 0xF8A733417D9161D9uL), (+1, -28, 0xEA916CB99B098ABEuL, 0x43F8E4F941FDDDA2uL)),
                ((-1, -32, 0x9CAE03112A9EA6F1uL, 0x830C0C4323A87CC7uL), (+1, -34, 0xA6D887800A097CFEuL, 0xB81AF9AFD9CA1FACuL)),
                ((-1, -39, 0x82CB554079B70D23uL, 0x27EF25AB19BC5B4DuL), (+1, -42, 0xF821910F7F209765uL, 0x93642D98BB1824E2uL)),
                ((-1, -48, 0xBBE1A5FD2DEFADFFuL, 0x929F38FAE7680D5AuL), (+1, -50, 0x99B2DC45839102ADuL, 0x7E9C87758E55E841uL)),
                ((-1, -58, 0x8B1C65BC465026B8uL, 0xABA6BF3F1A6D01B7uL), (+1, -61, 0xB38ADA78AB34123CuL, 0x16539866ED9B95B2uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0x94C8A9B79D8CC212uL, 0x1D974A4B5FA57A1BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -1, 0xD87175F1FE5831C0uL, 0x306E8E9EA3A370E8uL), (+1, -2, 0xAFF650330A25CD96uL, 0x753CCD556499E978uL)),
                ((-1, -3, 0x88CBD2C91D862756uL, 0x04011CF1C3C43C85uL), (+1, -5, 0xD2107AD67C2FD0A1uL, 0x790332C3C3DCC114uL)),
                ((-1, -7, 0xC53894CD418AACA3uL, 0x75EE4E97215B0EBCuL), (+1, -8, 0x8EDCE475DCC69F99uL, 0xF6F137A6F83A9152uL)),
                ((-1, -11, 0xB2F6867D3DC6760CuL, 0x8F3E130BB8FC888DuL), (+1, -13, 0xF4372A703A5B77EBuL, 0xA02A4413A9AF2AA0uL)),
                ((-1, -16, 0xD4E9A5E5A8BD1861uL, 0xCD402378F9823EC2uL), (+1, -17, 0x8883C2F21562FF73uL, 0x5A994174DB3842DCuL)),
                ((-1, -21, 0xA7C84FD46978F245uL, 0x5862514B8309D79DuL), (+1, -23, 0xC985C0F50D3FB791uL, 0x08CBCFD8FE1D7A37uL)),
                ((-1, -27, 0xAD187932FA329F9DuL, 0x3DEEBF04ABACD792uL), (+1, -29, 0xC1D4C4EEF4B6648FuL, 0xFA09F2AEC794FDE4uL)),
                ((-1, -34, 0xE1F3C422584640E4uL, 0xBDEDB24121D3C11BuL), (+1, -36, 0xEA54E9D09D3C6345uL, 0xE3725515989B8AD0uL)),
                ((-1, -41, 0xAF289901D75AE666uL, 0x38873BAF4CB8EDD9uL), (+1, -43, 0xA690F7B43FBD33E8uL, 0x8544872862F88356uL)),
                ((-1, -49, 0x904401F410732BF9uL, 0xE4F374DC7BE20117uL), (+1, -52, 0xF783CAE6BC360CDCuL, 0x62393805E00E7DD9uL)),
                ((-1, -59, 0xCC34D2DDB165B67BuL, 0x9FB08F4C3C396531uL), (+1, -61, 0x992BD8E2B6DB806FuL, 0x40FAC52101FC43F8uL)),
                ((-1, -70, 0x945BAFEF74D464EAuL, 0x4AADA7892CE0566DuL), (+1, -73, 0xB2BB6B2B462C50EEuL, 0x81CBD3759FFB0578uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xB4AF16B0E063702BuL, 0x4227083D9FF0D75CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -1, 0x81987178672C5080uL, 0x7571355150F97192uL), (+1, -3, 0xAFB407423E7431EBuL, 0xB86E01904BDD8138uL)),
                ((-1, -5, 0xA1A867A3DC106208uL, 0xB47CD67825F09AF2uL), (+1, -7, 0xD1709267300E03A6uL, 0x2A93739B5D398266uL)),
                ((-1, -10, 0xE626410F4D2A9723uL, 0x31363921CDC2FF52uL), (+1, -11, 0x8E37FE09B09DC47FuL, 0x535E7E24839F57C1uL)),
                ((-1, -15, 0xCE587B6FE7F393D8uL, 0xFFF2EE16D3657FAEuL), (+1, -17, 0xF2BB1B20C43BCF30uL, 0xDFE45F271B4AD1D9uL)),
                ((-1, -21, 0xF2A9E8FF99CBFEAFuL, 0x1AB430AC41A63A69uL), (+1, -22, 0x87772BF1FC7F9385uL, 0x0948F59B33E11FE6uL)),
                ((-1, -27, 0xBD171651D6DDFCA2uL, 0xACF7CAFEB7A0766EuL), (+1, -29, 0xC7A4879A1AF64801uL, 0x4CCBED0F1B3AE055uL)),
                ((-1, -34, 0xC0F079EC06217DC6uL, 0x124594A98594B5CDuL), (+1, -36, 0xBFB2B0CCEDC8E192uL, 0xD9AE6548DE5985D8uL)),
                ((-1, -42, 0xF91BC3F1CF93C10FuL, 0x610BCC3E89369F33uL), (+1, -44, 0xE75A37D319B06219uL, 0x4B0A23B3249BF043uL)),
                ((-1, -50, 0xBEF8830A6CB4651EuL, 0xC94357788FA7E9E3uL), (+1, -52, 0xA428D11432E03D46uL, 0x2E5164AB59606503uL)),
                ((-1, -59, 0x9B7614A7FDB0D805uL, 0xA94BEE465A3BC3D8uL), (+1, -62, 0xF380F3C200C285EEuL, 0x6129034616312EDAuL)),
                ((-1, -70, 0xD932B3D66A91F02DuL, 0xA0A4893C2E7CEF2DuL), (+1, -72, 0x966B22727CB89A12uL, 0x55ED4C167298F282uL)),
                ((-1, -82, 0x9B1E782529244850uL, 0x189B37D2F9EAF871uL), (+1, -85, 0xAF344EB6B2FBBD1FuL, 0xA053F32185A83C59uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xD3020CC1E31792BBuL, 0x4949593D383572D3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -2, 0x95E311D6BF8F8FEDuL, 0x12F84F52F921F524uL), (+1, -4, 0xAF54438376E8051EuL, 0x3492F5DED0253360uL)),
                ((-1, -7, 0xB935E3C929CCD23CuL, 0xA487D29F68A3A741uL), (+1, -9, 0xD08BE287EF973A7AuL, 0x9E378CD737C8127FuL)),
                ((-1, -12, 0x82A00158B2108E03uL, 0x8C7AC42E06E2953EuL), (+1, -14, 0x8D4EB9265767EA35uL, 0x2C60D1593B841732uL)),
                ((-1, -19, 0xE81ABA900DCD6948uL, 0xD85D687D0996E161uL), (+1, -21, 0xF0A7A9BD15D52CECuL, 0xAA7DD00F7524B75BuL)),
                ((-1, -25, 0x8740CC6CB4E545CFuL, 0xD71D6968CD7AF1E8uL), (+1, -27, 0x860428BB3D2E02F5uL, 0x80EB2629A5EEF7D2uL)),
                ((-1, -33, 0xD0E6EA614A7047D1uL, 0x7CCF98654477D5F8uL), (+1, -35, 0xC5143EB7EE62FBCFuL, 0x1C7D5015ECF71E58uL)),
                ((-1, -41, 0xD33DC603646D1AF1uL, 0x6DC65307B7402CFAuL), (+1, -43, 0xBCD3CAE6361B6D64uL, 0x480864EFB30AADB0uL)),
                ((-1, -49, 0x8720072D0E072B66uL, 0xB8527DB633D9ADF0uL), (+1, -52, 0xE365CA00C2529811uL, 0x853D5948DCE3FF92uL)),
                ((-1, -59, 0xCD38BCF64200D068uL, 0xA05BF3F557B2825EuL), (+1, -61, 0xA102555B8A6654C6uL, 0x2EFBC4DBB2333475uL)),
                ((-1, -69, 0xA560B8BBE6EEAF54uL, 0x974BEF282CFE25B0uL), (+1, -72, 0xEE53FE08AD90CE16uL, 0x794B371B2194E02DuL)),
                ((-1, -81, 0xE46925EACF20B19DuL, 0x994B54851026CD12uL), (+1, -83, 0x92EAD8781103D2F2uL, 0xAB34E877B3938866uL)),
                ((-1, -94, 0xA0A68E9358FD733AuL, 0x11BBEC5E314B8520uL), (+1, -97, 0xAAC82533A1CCE488uL, 0x6A50E63AA637B978uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm128_256 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xF06D47D6DBE46273uL, 0xC929EBC08EFCB1FAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -3, 0xA9930AE0646C5F79uL, 0xF1BC9AF55F9F2BFAuL), (+1, -5, 0xAEF7D6FB446AB5CFuL, 0x2A05AD75F44E4BE4uL)),
                ((-1, -9, 0xD00FDD38CD348DB3uL, 0x0F84E24F6DEDBF0AuL), (+1, -11, 0xCFB046E87A737009uL, 0x3D43B5E0E3DF6E95uL)),
                ((-1, -15, 0x91B61BAF0BAA0938uL, 0x9A75F49762033013uL), (+1, -17, 0x8C6FDEF03C811709uL, 0x9BC33D804B51197AuL)),
                ((-1, -22, 0x808BC29A48CBA2CDuL, 0x6CC7A5265FDB28E1uL), (+1, -25, 0xEEAEA6CC23519816uL, 0x08212E5D5C5B0CB7uL)),
                ((-1, -30, 0x94C18095E94D5024uL, 0x0972F04C860B2D89uL), (+1, -32, 0x84A58665FF0FE40FuL, 0xEE2DD0AB18A574ADuL)),
                ((-1, -39, 0xE41DD359E097499FuL, 0x1E02CF80766F90AFuL), (+1, -41, 0xC2AB7A0EBEC9B161uL, 0x4B5FCC91D48FF231uL)),
                ((-1, -48, 0xE4FDDBE2ACB68039uL, 0x0996E9247CCC5E3FuL), (+1, -50, 0xBA251574424FED15uL, 0x34CD65B2AA70EBCBuL)),
                ((-1, -57, 0x9161330344081D48uL, 0x54CF96DAE637F769uL), (+1, -60, 0xDFB939D6217073C3uL, 0x36A92FBE3586F413uL)),
                ((-1, -68, 0xDB0E075D91437E59uL, 0xE595D939DF8DAA96uL), (+1, -70, 0x9E197AA0F02328F9uL, 0x101995AC936E1F6FuL)),
                ((-1, -79, 0xAF04771A46589548uL, 0x88960066D9208D46uL), (+1, -82, 0xE99377902D326418uL, 0x529C7B86DD9F079AuL)),
                ((-1, -92, 0xEF5878C12D2233DBuL, 0x5AB5D8D6D6B3BFA6uL), (+1, -94, 0x8FB8BA4BD1C469B2uL, 0x015FCDB4E57D7ADDuL)),
                ((-1, -106, 0xA618EACED440290AuL, 0x7CF1B6ECAE2A58CCuL), (+1, -109, 0xA6C4B7686668E2D4uL, 0xA7A670CDFE2EB074uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm256_512 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 2, 0x86AA1257E3383901uL, 0xF189E2A88269F7BCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -4, 0xBCF198CBEB02A1E8uL, 0x3CBDB0D81638A64DuL), (+1, -6, 0xAEAB3E73C0DFED32uL, 0x40F584FE6EB07860uL)),
                ((-1, -11, 0xE694FD6DE641E4D0uL, 0x4B67269978AC3E02uL), (+1, -13, 0xCEFAEFBF4FC60642uL, 0xF5A5A97AC23D7607uL)),
                ((-1, -18, 0xA09BAB5637BED8FFuL, 0xA8B2E81863DDC755uL), (+1, -20, 0x8BB886420B1A5C6EuL, 0xB5D91BCDDFA05642uL)),
                ((-1, -26, 0x8CE94B49E22BCA7FuL, 0x1AD9F78824B60452uL), (+1, -29, 0xED10B33D63D8FDAEuL, 0x64EDD05695D35C64uL)),
                ((-1, -35, 0xA22835E2631542CEuL, 0x9F87711E6B9093A5uL), (+1, -37, 0x83872F2B72A4D04CuL, 0xCA9B607CB5EF0F8AuL)),
                ((-1, -45, 0xF740420FC4189044uL, 0xD0FF4A70C8DB9F7AuL), (+1, -47, 0xC0B5AF493CE91466uL, 0x1951917BB9992DDFuL)),
                ((-1, -55, 0xF6BDFE9E1989BD26uL, 0xA86FFD4D0564F964uL), (+1, -57, 0xB7F87B80DFBB367CuL, 0x55517FD868051A27uL)),
                ((-1, -65, 0x9BB065AEADA4F6DBuL, 0x36A196A9CC6D4930uL), (+1, -68, 0xDCC1B7A8FD5D0DA3uL, 0x04825DE2DD4DC596uL)),
                ((-1, -77, 0xE90F5A0D8AB3FC3EuL, 0xBD62EE7BF1A32D4CuL), (+1, -79, 0x9BC23DA089E3724CuL, 0xD0A9CE8B7E7AF3B3uL)),
                ((-1, -89, 0xB8DFC1740CA2C5F7uL, 0x87C6398122594F9CuL), (+1, -92, 0xE5C45C7149E78FCDuL, 0x81EA6526595BBB4EuL)),
                ((-1, -103, 0xFAB45A0B5A17EE4BuL, 0x8473024FA13FC806uL), (+1, -105, 0x8D2B5471F6B219A0uL, 0x05B89CA964B3C464uL)),
                ((-1, -118, 0xABF4EFB11AB32BAFuL, 0xDF710305B6FAA278uL), (+1, -121, 0xA393091C0D92F4F4uL, 0x23E2F72318489285uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm512_1024 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 2, 0x94F7FDD087C21280uL, 0xCBDFB96998EB06E1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -5, 0xD0283EB7593B3FDDuL, 0x943D732D0EA8F235uL), (+1, -7, 0xAE707A062A5B694AuL, 0xDFE856CD4DDF5FA6uL)),
                ((-1, -13, 0xFCF97A1B0AE3D477uL, 0x1F6E9A9B11A58D3EuL), (+1, -15, 0xCE702DEEB41FF829uL, 0xE9AF8C6FE7AB05C6uL)),
                ((-1, -21, 0xAF7615B5D30FA421uL, 0x072FD47CBDF39721uL), (+1, -23, 0x8B2C9C3987986798uL, 0x3D1A9030A3E22C7EuL)),
                ((-1, -30, 0x99479973F7541C5AuL, 0x6CC45DCC04DD9637uL), (+1, -33, 0xEBD5A872E16BF0BFuL, 0x53E0CE00A612C1C5uL)),
                ((-1, -40, 0xAF9C6EA9A2453E7EuL, 0x0B4F9EF671B71FD1uL), (+1, -42, 0x82ADD8A184EE1ADDuL, 0x874CF91082ECE649uL)),
                ((-1, -50, 0x8545AFF53AF592C4uL, 0x4FE3A1A9935B6B0CuL), (+1, -53, 0xBF39D4D79C172B3BuL, 0xFDA18022ACFA2349uL)),
                ((-1, -61, 0x845DDBBD1F5D91E7uL, 0xABF3A4D359F49EB9uL), (+1, -64, 0xB6544103032FC059uL, 0xD678FA7CE9D183F2uL)),
                ((-1, -73, 0xA634694E4126E5B2uL, 0x32DA0CC8CB1A12B5uL), (+1, -76, 0xDA85CBAAFD0E7AE6uL, 0x00CCE82F9EC5EBB7uL)),
                ((-1, -86, 0xF775FDFF7787AD9FuL, 0x2328D078BB3718A1uL), (+1, -88, 0x9A002A6AA751E954uL, 0x9234E9C935AC1FB4uL)),
                ((-1, -99, 0xC31E96F1D3050372uL, 0x4B7A2EC2ED628918uL), (+1, -102, 0xE2E9D2DD44B70694uL, 0x18023AA550329DF6uL)),
                ((-1, -113, 0x835A56AB78FAFCD9uL, 0x97E5EB58B921A588uL), (+1, -116, 0x8B4302C6DA128760uL, 0xA25523BEE5FBE64AuL)),
                ((-1, -130, 0xB25B9A34BBFC77F7uL, 0x62979612005288DCuL), (+1, -133, 0xA13162694E254DADuL, 0x2DC75CB081A97F78uL)),
            }));

            private static ddouble UpperValue(ddouble x) {
                Debug.Assert(x <= 0.5d);

                if (x >= 0.125d) {
                    ddouble y;
                    if (x <= 0.25d) {
                        y = ApproxUtil.Pade(x - 0.125d, pade_upper_0p125_0p25);
                    }
                    else if (x <= 0.375d) {
                        y = ApproxUtil.Pade(x - 0.25d, pade_upper_0p25_0p375);
                    }
                    else if (x <= 0.4375d) {
                        y = ApproxUtil.Pade(x - 0.375d, pade_upper_0p375_0p4375);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.4375d, pade_upper_0p4375_0p5);
                    }

                    return y;
                }
                else {
                    ddouble v;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 3)), pade_upper_expm3_4);
                    }
                    else if (exponent >= -8) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 4)), pade_upper_expm4_8);
                    }
                    else if (exponent >= -12) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 8)), pade_upper_expm8_12);
                    }
                    else if (exponent >= -16) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 12)), pade_upper_expm12_16);
                    }
                    else if (exponent >= -32) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 16)), pade_upper_expm16_32);
                    }
                    else if (exponent >= -64) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 32)), pade_upper_expm32_64);
                    }
                    else if (exponent >= -128) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 64)), pade_upper_expm64_128);
                    }
                    else {
                        v = Ldexp(RcpPi, 1);
                    }

                    ddouble y = v / x;

                    return y;
                }
            }

            private static ddouble LowerValue(ddouble x) {
                Debug.Assert(x <= 0.5d);

                if (x >= 0.125d) {
                    ddouble y;
                    if (x <= 0.25d) {
                        y = ApproxUtil.Pade(x - 0.125d, pade_lower_0p125_0p25);
                    }
                    else if (x <= 0.375d) {
                        y = ApproxUtil.Pade(x - 0.25d, pade_lower_0p25_0p375);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.375d, pade_lower_0p375_0p5);
                    }

                    return y;
                }
                else {
                    ddouble y;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 3)), pade_lower_expm3_4);
                    }
                    else if (exponent >= -8) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 4)), pade_lower_expm4_8);
                    }
                    else if (exponent >= -16) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 8)), pade_lower_expm8_16);
                    }
                    else if (exponent >= -32) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 16)), pade_lower_expm16_32);
                    }
                    else if (exponent >= -64) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 32)), pade_lower_expm32_64);
                    }
                    else if (exponent >= -128) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 64)), pade_lower_expm64_128);
                    }
                    else if (exponent >= -256) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 128)), pade_lower_expm128_256);
                    }
                    else if (exponent >= -512) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 256)), pade_lower_expm256_512);
                    }
                    else if (exponent >= -1024) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 512)), pade_lower_expm512_1024);
                    }
                    else {
                        return NegativeInfinity;
                    }

                    return y;
                }
            }

            public static ddouble Value(ddouble x, bool complementary) {
                if (x > 0.5) {
                    return Value(1d - x, !complementary);
                }

                return complementary ? UpperValue(x) : LowerValue(x);
            }
        }
    }
}