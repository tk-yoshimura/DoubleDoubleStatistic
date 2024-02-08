using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class WignerSemicircleDistribution : Distribution {

        public ddouble Radius { get; }

        private readonly ddouble pdf_norm, radius_sq;

        public WignerSemicircleDistribution(ddouble radius) {
            ValidateScale(radius);

            this.Radius = radius;

            this.radius_sq = radius * radius;
            this.pdf_norm = 1d / (PI * radius_sq);
        }

        public override ddouble PDF(ddouble x) {
            if (x < -Radius || x > Radius) {
                return 0d;
            }

            ddouble pdf = 2 * this.pdf_norm * Sqrt(radius_sq - x * x);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (x <= -Radius) {
                return (interval == Interval.Lower) ? 0d : 1d;
            }
            if (x >= Radius) {
                return (interval == Interval.Lower) ? 1d : 0d;
            }

            ddouble v = x * this.pdf_norm * Sqrt(radius_sq - x * x) + Asin(x / Radius) * RcpPI;

            ddouble cdf = (interval == Interval.Lower) ? 0.5d + v : 0.5d - v;

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (Abs(p - 0.5d) < 1e-31d) {
                return 0d;
            }

            if (p > 0.5d) {
                return -Quantile(1d - p, interval);
            }

            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble u = Sqrt(2 * p);

            ddouble quantile = QuantilePade.Value(u);
            quantile = (interval != Interval.Lower) ? quantile : -quantile;
            quantile *= Radius;

            return quantile;
        }

        public override (ddouble min, ddouble max) Support => (-Radius, Radius);

        public override ddouble Mean => 0d;
        public override ddouble Median => 0d;

        public override ddouble Variance => radius_sq / 4;
        public override ddouble Skewness => 0d;
        public override ddouble Kurtosis => -1;

        public override ddouble Entropy => Log(PI * Radius) - 0.5d;

        public override string ToString() {
            return $"{typeof(WignerSemicircleDistribution).Name}[radius={Radius}]";
        }

        private static class QuantilePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0_0p0009765625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0x8000000000000000uL, 0x0000000000000000uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 22, 0x939DC33169785A61uL, 0xE18F2773FEAD153EuL), (+1, 22, 0x939DC332F7F5B676uL, 0x93985F15E1F5B3CBuL)),
                ((+1, 42, 0xA95919784CD0AF37uL, 0x95F46DB68ED24D5DuL), (+1, 42, 0xA9591987590DD79BuL, 0xFFCC96E13BD0F728uL)),
                ((+1, 61, 0x9E03E6A9B7390AE0uL, 0x3A5EFDA53CE6A114uL), (+1, 61, 0x9E03E6DC54FBC1A9uL, 0xDB05F554D1301F52uL)),
                ((+1, 79, 0x944C223FF57B06E1uL, 0xF7C8A21BDD72CCB9uL), (+1, 79, 0x944C22BB908D0C45uL, 0xE9C13C2D91167B52uL)),
                ((+1, 96, 0x9DA71D1EEDC4338BuL, 0x90722D649630D9A5uL), (+1, 96, 0x9DA71E3DF7095164uL, 0x5BCC52A6807F1164uL)),
                ((+1, 112, 0xCCC00CE0D6591254uL, 0x8513D68DA6C8CF1DuL), (+1, 112, 0xCCC00FBB22CBF27EuL, 0x2AD15FB17785836BuL)),
                ((+1, 128, 0xAB10E78FB98C8FF8uL, 0xC20221E3DFC2BC57uL), (+1, 128, 0xAB10EBE620D6A718uL, 0xDB94C7B201DA3FB5uL)),
                ((+1, 143, 0xBEDA2654698212F0uL, 0x1BC79DC163B387A2uL), (+1, 143, 0xBEDA2EAFA6FB1B9BuL, 0xB50595D5ADBB025BuL)),
                ((+1, 158, 0x921EF2E5AE22909FuL, 0x90A8CF3E984C4B29uL), (+1, 158, 0x921EFD86F060F73EuL, 0x674E76A63B2AFFABuL)),
                ((+1, 172, 0x9CBA9BAA292B9955uL, 0x57AAAA31D8B11FC9uL), (+1, 172, 0x9CBAAE121EA37545uL, 0x9DDAD90F414D1EF5uL)),
                ((+1, 185, 0xEF23BC7E966DFAF3uL, 0x16BCFF50AABFAE16uL), (+1, 185, 0xEF23E8DE46309A5DuL, 0x3B2900C6D291638CuL)),
                ((+1, 199, 0x8338925F1574E05FuL, 0x4C3D71A6FA3AAA7AuL), (+1, 199, 0x8338B83DEDD3E8C6uL, 0xC3F3A8FE48BD8457uL)),
                ((+1, 211, 0xD0C09D6B084F97A7uL, 0xEF63A092659D19CBuL), (+1, 211, 0xD0C0FA18C25FFB78uL, 0x4853369FFCADC569uL)),
                ((+1, 223, 0xF1D0433E71A41F28uL, 0xBEDB47A85D3E858FuL), (+1, 223, 0xF1D0E7458B733973uL, 0xDC1763FC9C83190BuL)),
                ((+1, 235, 0xCC51B6C5E973FF6DuL, 0x5C3E52B948835647uL), (+1, 235, 0xCC5289E71EB7051EuL, 0x2D3992E62B77EE3DuL)),
                ((+1, 246, 0xFB93CA434BA38066uL, 0xC75B1B613C35415CuL), (+1, 246, 0xFB95569C02526023uL, 0x7AB2A6D035F4B950uL)),
                ((+1, 257, 0xE0C8C721376271EDuL, 0x6C5854D5540290BCuL), (+1, 257, 0xE0CAE5A27538114FuL, 0xEAAC7DF5C61E7F11uL)),
                ((+1, 268, 0x90ABBF6783167131uL, 0x9994232E023FC4D0uL), (+1, 268, 0x90ADDB2C4E2C9C76uL, 0x4900BBAB02FE70DDuL)),
                ((+1, 278, 0x84A180CD3474DDD8uL, 0xF504BA4395160C47uL), (+1, 278, 0x84A488B5095049EAuL, 0x9821648183C8A2B7uL)),
                ((+1, 287, 0xAA6F641A6DF41213uL, 0xDFB3D0E20266206CuL), (+1, 287, 0xAA759FA35FD121CAuL, 0xE6D76ABA1E5302F9uL)),
                ((+1, 296, 0x9617BC313DA23AB2uL, 0xE0F4D00DCDC20BDAuL), (+1, 296, 0x9620C57DC94A4FE6uL, 0x549038781589552EuL)),
                ((+1, 304, 0xAF9CF6E9AD58FA33uL, 0xA26A941F55C9369BuL), (+1, 304, 0xAFAF125AD4683053uL, 0xF5209F8F7B5AEB2CuL)),
                ((+1, 312, 0x82A047E9ED94ADF2uL, 0x275223D7D8D19DD2uL), (+1, 312, 0x82B8AA77D91A65DBuL, 0xA3D8454823BC96FEuL)),
                ((+1, 318, 0xE7ACF9B4E71A020DuL, 0xDFEE34809C8A3D53uL), (+1, 318, 0xE801E130F2017478uL, 0xE40E12324C4F9DB4uL)),
                ((+1, 324, 0xDD85AA3A7E23BDEBuL, 0xD18615AAC5CD9FC8uL), (+1, 324, 0xDE3A4CA987449C56uL, 0xD0A3792592E41B12uL)),
                ((+1, 329, 0xBF9E278677E2FA22uL, 0x29B3CDE0448DDF7EuL), (+1, 329, 0xC14D647196F62F23uL, 0x0141B4A40DB30EA2uL)),
                ((+1, 332, 0xCB13A6F467DB6CC7uL, 0xFAAB3E3CE5BB67A2uL), (+1, 332, 0xD2E9858C93322E45uL, 0xA1989C2BFFB03E87uL)),
                (0, (+1, 331, 0xD78B2EE2BAA62B48uL, 0xE598EC2D23640790uL)),
                (0, (+1, 332, 0xF357FA18EE663888uL, 0xA84C7CA7BCA79E1FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p0009765625_0p001953125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFFFA60DDCD90C4BFuL, 0x30D01536C82D0EB1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 12, 0x83731AC7CD8C8F1FuL, 0x32E3BBDAF1263C10uL), (+1, 12, 0x8376EDAE6E8589B2uL, 0xC06C380DBC373990uL)),
                ((+1, 22, 0xE1D5E4F22337B5C2uL, 0x6F97DD0F59107467uL), (+1, 22, 0xE1DEDC14E904EDE3uL, 0x135F8246F5E65598uL)),
                ((+1, 32, 0xD2E6B74BE7941870uL, 0x0D7D54A461391DE0uL), (+1, 32, 0xD2F2918C71026F4AuL, 0x0A6AFD434858F575uL)),
                ((+1, 41, 0xE91CC10FA2CAEA3FuL, 0x1A7694A25F88E876uL), (+1, 41, 0xE9302D16F0EF2BF9uL, 0x399915A02859BB50uL)),
                ((+1, 50, 0x9BB79EE2B00B6A20uL, 0x4B1BDA4E64227E9AuL), (+1, 50, 0x9BCC0AC90B6482B5uL, 0xF757BD15F2B09FA9uL)),
                ((+1, 57, 0xF6923C6C7ABA447DuL, 0xE22542F183D0AC76uL), (+1, 57, 0xF6C970B0E6B61EFFuL, 0xE3FD33577041F8D5uL)),
                ((+1, 64, 0xD9AF5C8EA0362AEAuL, 0xBCFAEEF3D3F3FE35uL), (+1, 64, 0xDA0CDECDF943980DuL, 0xD968E7E7F4F2089AuL)),
                ((+1, 70, 0xBD462EF83DCD5744uL, 0x84FB819D6086CFE2uL), (+1, 70, 0xBE01BCF6C0736A89uL, 0xB3456ACEC8226AFFuL)),
                ((+1, 74, 0xF9301A03EC13EB01uL, 0x30AFC3668268748BuL), (+1, 74, 0xFC521EF6D20979B8uL, 0x630A900320F9D1FDuL)),
                ((+1, 76, 0xDC04A23731B0368FuL, 0xEB680BE53F774315uL), (+1, 76, 0xF3560590D4F2C6B2uL, 0x9900A78439E262CAuL)),
                ((-1, 76, 0xB94D9F066D374226uL, 0xE336A67324C667A8uL), 0),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p001953125_0p00390625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFFF1D590367F0AF5uL, 0xC83D456191F9A48CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 11, 0x8A39BEB9AC4E14C6uL, 0x11F3B89EB7510D88uL), (+1, 11, 0x8A43C1C5124475B3uL, 0xEA31C06E850551F4uL)),
                ((+1, 20, 0xFBE95B6FF7942921uL, 0x5AD68783AE2AA6DEuL), (+1, 20, 0xFC01E5A0764796E1uL, 0x3C566C9B26EDDEA0uL)),
                ((+1, 29, 0xFC50315E33E9EC54uL, 0xD2CF738A4D0494B8uL), (+1, 29, 0xFC725E5DA69A794EuL, 0xC5013EE556544A31uL)),
                ((+1, 38, 0x97B8C46EDF579F37uL, 0x147CE98DED295665uL), (+1, 38, 0x97D68DA9F29A5E69uL, 0xFB3C4C86D02D9882uL)),
                ((+1, 45, 0xE0DBF05BC1A21BB0uL, 0xE0C4968453EE9B68uL), (+1, 45, 0xE11F6CCE5B09E51CuL, 0xD2C33476E51AC3B4uL)),
                ((+1, 52, 0xCAEA75838E65AE80uL, 0xB74A241740FB1277uL), (+1, 52, 0xCB4E6E445833D7A2uL, 0x8EEE56C48E4490F1uL)),
                ((+1, 58, 0xD486EDC30A3CAC3EuL, 0x78CFE1C698B9B44CuL), (+1, 58, 0xD544FE5CBA8F29A4uL, 0x80EA79A083FD19F4uL)),
                ((+1, 63, 0xE9BCD60CDFBCBBADuL, 0x614BE9D57969D511uL), (+1, 63, 0xEB77D97E46831496uL, 0xE24CC2E43BB078ECuL)),
                ((+1, 67, 0xDAF1867408E8DFF5uL, 0xA9E32D48AA5F13C8uL), (+1, 67, 0xDF7BEC5A4D4B529CuL, 0xDDA3C41BE52A4CAAuL)),
                ((+1, 69, 0xB96F93F5D94F3919uL, 0xE45EDC88AB17D5F1uL), (+1, 69, 0xD04298F23F065EF0uL, 0x133CC85A48D89D22uL)),
                ((-1, 70, 0xAFD35EB2F4F46EAFuL, 0x7BADF042AD6A8FA5uL), (-1, 69, 0xC86CE9B8C6351572uL, 0x019203DD04A4E950uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p00390625_0p0078125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFFDC4DA6E63FB4EBuL, 0xA683B0F635283B32uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 10, 0x8A308187AE7D909CuL, 0xD5089B6F1C6F5512uL), (+1, 10, 0x8A49BD1807819047uL, 0x754940001F447AD0uL)),
                ((+1, 18, 0xFBDBC2A574B1E179uL, 0x542C348B392CA881uL), (+1, 18, 0xFC199E00C34E8DFBuL, 0xAB9ABD1637129DF8uL)),
                ((+1, 26, 0xFC433BE825A6FD38uL, 0xC75FC53468A675BBuL), (+1, 26, 0xFC9965804F1101FEuL, 0x405E6A1EA912BBF9uL)),
                ((+1, 34, 0x97ADF53019F4C716uL, 0x22C2AC24D780ADFAuL), (+1, 34, 0x97F912085DA65DFEuL, 0x7B9DC13D5CC67A1FuL)),
                ((+1, 40, 0xE0BC41002752E116uL, 0x1C287BC095DEDF72uL), (+1, 40, 0xE1667C1CD8D33D28uL, 0xEC59508689292796uL)),
                ((+1, 46, 0xCAA82DA497A472B9uL, 0x372F50E379A081B0uL), (+1, 46, 0xCBA47183A3A7F5F9uL, 0xD808FF1D9FDD0106uL)),
                ((+1, 51, 0xD3DA6DA6242A0DF2uL, 0x4320F28528B648E2uL), (+1, 51, 0xD5BA419092427DC3uL, 0xE1CE5BC3738F4A03uL)),
                ((+1, 55, 0xE7B7486EE0D34E3CuL, 0xE74C00CC6381A96FuL), (+1, 55, 0xEC1666F4E930FACBuL, 0xE43E14D49E0AD8C2uL)),
                ((+1, 58, 0xD488DB87F77196EDuL, 0x39A53CFE56DB11E2uL), (+1, 58, 0xE003ED816BA59847uL, 0xF10A983038BDE196uL)),
                ((+1, 59, 0x94C98D1B78FF3815uL, 0x5203B8F25E0B0143uL), (+1, 59, 0xCE945C76086FCAF3uL, 0x365634D7EBBD5CB2uL)),
                ((-1, 60, 0x98FAB946D3EB906AuL, 0x370ABDEBA6826B02uL), (-1, 58, 0xE407EA7A8D07BFC2uL, 0xF9176E50D1A1142DuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p0078125_0p015625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFFA60AE60D349F20uL, 0xDA3BAE372B0B8312uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 9, 0x8A318935BDBE7FF8uL, 0xD92FA4EADEAB67DEuL), (+1, 9, 0x8A712DD0F606C851uL, 0xE2400B2A9D8A9AC6uL)),
                ((+1, 16, 0xFC1BABE471BE6F01uL, 0xFE9968B8A251DF29uL), (+1, 16, 0xFCB7D9BB419BBA00uL, 0xE42B29584F70F5ACuL)),
                ((+1, 23, 0xFCC7A9784F60B350uL, 0x673F1C803BAF0B5AuL), (+1, 23, 0xFDA17A746A536E2BuL, 0x77B5F88C610D694CuL)),
                ((+1, 30, 0x9828B1C3C984C294uL, 0x445FEF1F51FE2C8DuL), (+1, 30, 0x98E6DD3E5375BE3EuL, 0xCC5E143F049023EAuL)),
                ((+1, 35, 0xE1AC016F28E72667uL, 0x8A34AAC0B2AFC487uL), (+1, 35, 0xE35BCEA0D7BA29A4uL, 0x8E82CF013370B480uL)),
                ((+1, 40, 0xCB96BAA28DCD17A2uL, 0xF82F22B7166DFB0AuL), (+1, 40, 0xCE1827B044A3A360uL, 0x8F82C003CA7E28D0uL)),
                ((+1, 44, 0xD479F77C2233DC40uL, 0x4C94199C89E5800CuL), (+1, 44, 0xD941E00130B9FDD8uL, 0xE1656A9E45053F5AuL)),
                ((+1, 47, 0xE63656A9ACDAA6F6uL, 0x3AB0FAD317F1522CuL), (+1, 47, 0xF1694F1C3FF4D654uL, 0xD6BEDB409C16A2A3uL)),
                ((+1, 49, 0xC915D3266F6FE984uL, 0x48DAB32DA259E109uL), (+1, 49, 0xE6AE93803DC49E43uL, 0x2A989B1C2BB3BC4FuL)),
                ((+1, 47, 0xEEA54E587B525FADuL, 0xC7AC9939F778C46CuL), (+1, 49, 0xD21E92024A443E0AuL, 0x89DC9B388B6029A1uL)),
                ((-1, 50, 0xA75F62A9B1E04E73uL, 0x0A930AE6553F621DuL), (-1, 48, 0xA0F069D43BAB47E8uL, 0x5CC8FEF0ACDF0193uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p015625_0p03125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFF1D462DFD261EF3uL, 0xD83A8C47FD7E32A2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 8, 0x8C148759BF007292uL, 0xA2E50D284E2E961BuL), (+1, 8, 0x8CB6F3FDCF900D83uL, 0x4CEFB79E0B55A2D5uL)),
                ((+1, 15, 0x822268CE48D27DCCuL, 0xACCA57658C4580F5uL), (+1, 15, 0x82ECA1DD29804AE5uL, 0x062B2CF6C1206407uL)),
                ((+1, 21, 0x85541A0E8F3E766BuL, 0xA7B2A916E58FACEAuL), (+1, 21, 0x867306F2C22BD280uL, 0x1C598504DB7459ECuL)),
                ((+1, 26, 0xA4B09F3B220A87FFuL, 0x80B72C065F3F6665uL), (+1, 26, 0xA6B009E35A81BE48uL, 0xEEB9E9C70ECADD8CuL)),
                ((+1, 30, 0xFBE45D85405C7823uL, 0x5DB5B873BA5328F2uL), (+1, 31, 0x804573F618E5A6B4uL, 0x2B2DF2F6FF75671CuL)),
                ((+1, 34, 0xEBC248DAECD65FFAuL, 0x2AEBFF9DFE4B6FF1uL), (+1, 34, 0xF2E21C2C0D3A8E26uL, 0xA84B809D8A28B0C3uL)),
                ((+1, 38, 0x805EC293F9CBA29CuL, 0x709DE231D14A068DuL), (+1, 38, 0x876FD2576E4A076AuL, 0x8D478925DA00E83AuL)),
                ((+1, 40, 0x90D45FB49E709F0DuL, 0x8B459E6D92FCAECCuL), (+1, 40, 0xA23DAE3CCCD9B34CuL, 0xAC6FB26D65B61B68uL)),
                ((+1, 40, 0xF508D19EF2D2F49FuL, 0x8FAA05B766C416EEuL), (+1, 41, 0xABC545B2D75EC706uL, 0xA49D741697D6010DuL)),
                ((-1, 39, 0xC9FFA03F5D631D42uL, 0xE106E22B648FACB8uL), (+1, 40, 0xAEAEECF0D24E6A5EuL, 0xE4688D14EA09C1C4uL)),
                ((-1, 41, 0xA7549FA03AF7DBEBuL, 0x7AD527EED23E0BDCuL), (-1, 39, 0x81AEFD58C5B33447uL, 0x413531E8AFFE042BuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p03125_0p0625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFDC462A25AD22C1AuL, 0x525AAFD7BD7C11C0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 7, 0x8C55F030EB38AEA2uL, 0x456E9B5C3B5016D3uL), (+1, 7, 0x8DF23CFF8EB485A7uL, 0x045DF985C17FB873uL)),
                ((+1, 13, 0x832DBF450D64D37BuL, 0x0D2545033D8F687EuL), (+1, 13, 0x8532A54DBF5CA74FuL, 0x8E54F3BAFDDF9B12uL)),
                ((+1, 18, 0x86F2E65B082012D9uL, 0x399D0C8B333F21FBuL), (+1, 18, 0x89D54F3CA19F3076uL, 0xFC5475AE05255C01uL)),
                ((+1, 22, 0xA6A92589ACA007E9uL, 0x7E38C65BE70B8B74uL), (+1, 22, 0xABD5724B9587C624uL, 0xA44F44E9F24AFE93uL)),
                ((+1, 25, 0xFCA299FEB5441BCFuL, 0x7DDF0B39976ADA57uL), (+1, 26, 0x845DE844D7B3E977uL, 0x18EAACC0AE3BE7BCuL)),
                ((+1, 28, 0xE5EC564E8A06F09BuL, 0x34145701EA9302F2uL), (+1, 28, 0xF87B602A52339DDFuL, 0xE00943F39F0E924CuL)),
                ((+1, 30, 0xE82E53BD41BDE9F8uL, 0x99200CEECC1DC0FDuL), (+1, 31, 0x86656EC313E993D1uL, 0xD34B7489155B58F2uL)),
                ((+1, 31, 0xCE756C877D3A055CuL, 0x6D48A30918768E26uL), (+1, 32, 0x9369323E32FB383DuL, 0x83428F23B143E361uL)),
                ((-1, 27, 0x83189F006AAB4AA6uL, 0xF784953377B46C09uL), (+1, 31, 0xE3EC85127BF789DDuL, 0x509A808B12A48250uL)),
                ((-1, 32, 0xA602FA871BCC7B17uL, 0x3C2815705DF35CA7uL), (-1, 29, 0xD4A5608470B158E0uL, 0xF42220FC7923156BuL)),
                ((-1, 30, 0xBB055EA211419F73uL, 0xD19F2DFA684D74C6uL), (-1, 29, 0xAFAB105A41B54DEBuL, 0xF4CFAB7FD1B071B5uL)),
                ((+1, 29, 0xC879749CE4160E80uL, 0x4B465628532CC599uL), 0),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p0625_0p125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xFA5DB25113B8F97CuL, 0x1A2B2841C9D15774uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 6, 0x9BEFE6A35705EA7DuL, 0x38FF789745DEE210uL), (+1, 6, 0xA0688F483136F0E8uL, 0x2A50ACEB7572BD8CuL)),
                ((+1, 11, 0xA5B346B7F83B006EuL, 0x72F5E0D28862961FuL), (+1, 11, 0xABEC2794062296E1uL, 0x7574A4BBAE6CFA7AuL)),
                ((+1, 15, 0xC3C3CF0BD9BD220BuL, 0x24C8DE1D3F8C102FuL), (+1, 15, 0xCDBB5973D5896FA1uL, 0xC233032C072F0443uL)),
                ((+1, 19, 0x8C4D509F40657198uL, 0x546F101A6A138EADuL), (+1, 19, 0x967179D721614DF4uL, 0xC37BE8724F3CC6A6uL)),
                ((+1, 21, 0xF91863CA993CFFC4uL, 0xD2EE953D8DF71F94uL), (+1, 22, 0x8A30D6F1AD2FC9D8uL, 0x84B57C439AB862B5uL)),
                ((+1, 24, 0x84E4FFA8734F1550uL, 0xA4069818236BFD17uL), (+1, 24, 0x9D5BD81888016F9AuL, 0x07F777A230BF0C7FuL)),
                ((+1, 25, 0x97B12D6FC96757D1uL, 0xB00B957E2380190CuL), (+1, 25, 0xD11DC94898EDA62AuL, 0x75CBA8B247143B48uL)),
                ((+1, 24, 0xD67FB476D6C3DB99uL, 0x39FB94A09BB1A6CDuL), (+1, 26, 0x89757BBDDD89A8BFuL, 0xC9D259705F9D42FEuL)),
                ((-1, 25, 0xC3977A1A72201A9AuL, 0xB922509B3078A35DuL), (+1, 24, 0x96CDEA1E004CB346uL, 0x23303651766B8D09uL)),
                ((-1, 26, 0x95B41F2081433DF1uL, 0xE6136540C9FE80DAuL), (-1, 25, 0x839213D736F0464BuL, 0x1DE1164653FD74CDuL)),
                ((+1, 21, 0xF5C5C419C6C44A1EuL, 0x9D9BF2F639B3C16EuL), (-1, 23, 0xC5E1AE9FFB556811uL, 0xD14BE8FD8B8619EEuL)),
                ((+1, 24, 0xAE206391D8D2A3C0uL, 0xBF921B84A95B7925uL), (+1, 21, 0x9D823A8E0B0FEBB5uL, 0xB39C595E9B49F0EBuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xF1C137EE8FD17E41uL, 0x43225FDE4CDB0AE5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 5, 0x80ABFA566830D29DuL, 0xE43549FE6584E51AuL), (+1, 5, 0x8AC8396F81F10FE0uL, 0x0966102445D671FCuL)),
                ((+1, 8, 0xE24C99DD2A764AB8uL, 0xB4DA07F1061FA863uL), (+1, 8, 0xFB0A5112E543D0D3uL, 0xA881439090B5EB41uL)),
                ((+1, 11, 0xD262FDF79ADCDBF2uL, 0x098D001ADFB7907FuL), (+1, 11, 0xF471CA109E58F88AuL, 0xA5ECFC92B08DF695uL)),
                ((+1, 13, 0xD8E15BCB0F2370B0uL, 0xB68B0EDC0852ED8BuL), (+1, 14, 0x8943EB87AB3EB8B3uL, 0xCA85232B3AC2E3EEuL)),
                ((+1, 14, 0xE18D4AE19B79A350uL, 0x2FFF1A8E6190EF4DuL), (+1, 15, 0xAE23803C1BB6C36CuL, 0xC72CDF1083C46487uL)),
                ((+1, 13, 0xD047F3A6F2781B42uL, 0xAEE65602D9A017CEuL), (+1, 15, 0xD3BB2E01487520E6uL, 0xB62BBBF551CB307CuL)),
                ((-1, 15, 0xCFABD9E9B8C40D02uL, 0x33249401C6837332uL), (+1, 12, 0x900226D755D29B02uL, 0x651837A1CFC47A31uL)),
                ((-1, 16, 0x8F18FA0D0B504EBEuL, 0x3FA13A57D024324FuL), (-1, 15, 0xBE3C8E8758D55E3CuL, 0x3C3977A76F6B734FuL)),
                ((+1, 13, 0xB109FD9541B90A02uL, 0xF2373032002CECA6uL), (-1, 14, 0xB54315433676E2B2uL, 0x8BFF229DAA9356E9uL)),
                ((+1, 15, 0xB2BCC91441D77E82uL, 0x40C16E3A481E17DDuL), (+1, 13, 0xB1378ACCD4DBA5C3uL, 0x9AD779D2651C3632uL)),
                ((+1, 10, 0xC6D7DC2B87685942uL, 0xEF6E4BFB327C0F07uL), (+1, 11, 0xF8A951916449DDF3uL, 0x7059889D324C0BC8uL)),
                ((-1, 12, 0xCC2201C537F3A76BuL, 0x9D672410E20D6379uL), (-1, 9, 0x875FBC2551504B62uL, 0xF2C694AB8DAFC3FFuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xDBC8CEEC425E1B00uL, 0x8C18E50AB00815FFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0xA07D095D798E1601uL, 0x684E19AB309608F5uL), (+1, 3, 0xC934C85CF2AF07A2uL, 0xCCE7ACA1110F5821uL)),
                ((+1, 5, 0xA10931FCF3D76468uL, 0xDFDC5402E4691A54uL), (+1, 5, 0xEB1E0CCE7104036DuL, 0xBC8F57294D17DF48uL)),
                ((+1, 5, 0xC2946927638C1BF5uL, 0x39E4FF4D58919969uL), (+1, 6, 0xEA2F9858AFA44B31uL, 0xE56D587DBA922393uL)),
                ((-1, 6, 0x904F4DCF229CB076uL, 0x5C54A2620F679FD4uL), (+1, 5, 0xDDAD409FFBA0C582uL, 0x322E58CD087980EEuL)),
                ((-1, 7, 0xC44A3845FAB20893uL, 0xB7BCC9D762A94C9EuL), (-1, 6, 0xE7BBEAE426570674uL, 0x087025AF256E7507uL)),
                ((-1, 4, 0xDA749097D2203BA4uL, 0xB37DCAF2C979A666uL), (-1, 6, 0xE1E39BC23E8984AEuL, 0x6602F6A918FC69E5uL)),
                ((+1, 7, 0xA2DF4CADBC2E7D2CuL, 0x55C50F8B7ADF763BuL), (+1, 4, 0xF068BF31508108BEuL, 0x79B169B9C23899A1uL)),
                ((+1, 5, 0xB2BDBE378DCCD46AuL, 0xC3E0818BE1C04531uL), (+1, 5, 0x91B71CFFB4C2B435uL, 0xACA8564CF89EF3FAuL)),
                ((-1, 5, 0xAB2AC1325F253786uL, 0x5AC63B73141A8CF9uL), (-1, 1, 0xC9202BE28BAF8FFAuL, 0xBDD98BD223019949uL)),
                ((-1, 2, 0xCA3777817CF9D5DDuL, 0x94EC88C4D994F9D5uL), (-1, 0, 0xBC7A826F7CE69EF1uL, 0xC6545E256A83BB5FuL)),
                ((+1, 0, 0xE6AA0D47629D1E5FuL, 0x640148054FC1DF9DuL), 0),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xC11E903B82FC24F3uL, 0xF7922281D68BA18EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0x8808AD99C847A2E2uL, 0xE59D966DB07B85AEuL), (+1, 2, 0xDA642DFC75C17BFEuL, 0x01C7D771F581FDF1uL)),
                ((+1, 1, 0xFC17438453CC7A0DuL, 0x9D92E23A5A51FF1BuL), (+1, 3, 0xDFDC0D5EEC9321FDuL, 0x0DE32A815945C196uL)),
                ((-1, 3, 0xF056D03E5422868BuL, 0xC4756F6E85132507uL), (+1, 0, 0x81D12A84A214DA74uL, 0x76B8666F5E93B252uL)),
                ((-1, 4, 0xCD17BA2B534C91BEuL, 0x91FC4D213251A6B4uL), (-1, 4, 0xC336A0316EE78150uL, 0x4C2C96CD430A0319uL)),
                ((+1, 3, 0xB2F3C7D9868B24BFuL, 0xDF927129D90CEAF0uL), (-1, 3, 0xDADF137A46251F35uL, 0x9C4EA6C0F987E705uL)),
                ((+1, 5, 0x8073C4B1A82E634BuL, 0xE85845DF0281095EuL), (+1, 3, 0xD6A76716A8E434DEuL, 0xDAB0651954E06D35uL)),
                ((-1, 0, 0xAC35FEA27B0A7966uL, 0xD560E1BD6D50F276uL), (+1, 2, 0xF816144A36045AE6uL, 0xF99DE0B2686BB8FAuL)),
                ((-1, 3, 0xD622F44974483A1BuL, 0xE4126BCC1CD9691AuL), (-1, 1, 0xAFA7BF9BEDFFDEA4uL, 0xFEF562399D31B1ADuL)),
                ((-1, -3, 0xC6BFE8F3C0EFE2E0uL, 0xC94414FAF52FC397uL), (-1, -1, 0xDDB227936D09B080uL, 0x33172CEB37ED945BuL)),
                ((+1, 0, 0xB73AE0C510064179uL, 0x2F97A54D55178AE8uL), (+1, -4, 0xD5F4494A42738994uL, 0xB03E3758BA55530EuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p5_0p625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA27C0011BFAB4AFCuL, 0x8BD50B5B686256EBuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xF04C90EA94D43982uL, 0xAA4DF71D3BD788F4uL), (+1, 1, 0x95CF1FDA173BA7E9uL, 0x6667A6CA0579D843uL)),
                ((-1, 2, 0x91D24E0A49F4C8ACuL, 0x3A1659521B93EB24uL), (-1, 1, 0xAC7970528D829F5CuL, 0xE7231B90F1850446uL)),
                ((-1, 1, 0xBF597EF0CFE73293uL, 0xA51E9BEB992A2F60uL), (-1, 2, 0xEA63A9265E46BFB1uL, 0x7B09CCC952EBBF0AuL)),
                ((+1, 3, 0xA7D03C6AF186E153uL, 0x8E72CDA698F62802uL), (+1, 1, 0xC29417B9A59C1FAEuL, 0x70127AE411566649uL)),
                ((+1, 2, 0x9B0D42278A445123uL, 0xD2392B705959F58FuL), (+1, 2, 0xF6FA813678A806DAuL, 0xCEAC8C8C431718A4uL)),
                ((-1, 3, 0xA31F5710BCF95D0BuL, 0x503547F00276FE44uL), (-1, 1, 0x816B7C944479F649uL, 0x1C26EBBC7010AA75uL)),
                ((-1, 1, 0xA0875DB020014036uL, 0x046A12F85C0F72C9uL), (-1, 1, 0xBA58AEB203CD29CEuL, 0x907DD0E6293A34EFuL)),
                ((+1, 1, 0xFB605C145C4595B5uL, 0x2F6120081AB8093AuL), (+1, -1, 0x9BD1FCD856494300uL, 0xB5FB755E9324DA48uL)),
                ((+1, -2, 0x9B882D87B96A87FEuL, 0x105CC4A32EEBDA80uL), (+1, -2, 0x860700FAC2B59547uL, 0x3E587E9399DC9AABuL)),
                ((-1, -2, 0xCD43F27261D83465uL, 0x373DC6CEA66E30A8uL), (-1, -6, 0xE9F4891C52A5649DuL, 0xB7F9CEA6D892BF27uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p625_0p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x80166C22E1C6A174uL, 0xDAE99E6883AFFC61uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -1, 0xFAF03C8B64E14549uL, 0xE4F2B473D39731DCuL), (+1, -2, 0x9D3C3A3A5221AAA5uL, 0x24D7BA0DA4A55B62uL)),
                ((-1, 1, 0xACD9E07014CB5EE2uL, 0x0D545CBD2558B5A1uL), (-1, 1, 0xEFC68CFF96984193uL, 0xDBF0F33D326711FCuL)),
                ((+1, 1, 0xF077227E76018BE4uL, 0x80DE93081E566367uL), (-1, -1, 0x8BE8A808FD249552uL, 0xA60BDD74AB9E8D70uL)),
                ((+1, 2, 0x90790BA693B6C593uL, 0xF7AEB4976E1F5D9DuL), (+1, 2, 0x92F0AF62DF6C20E4uL, 0x6592374339C9D9E9uL)),
                ((-1, 2, 0x974BC6F9C508522FuL, 0xC749E801CCFB1535uL), (+1, -4, 0xC2FF88DAB7293D9DuL, 0xB45507AB00008EADuL)),
                ((-1, 1, 0xB136E6F9C2152995uL, 0x7BDBADD91FF7F895uL), (-1, 0, 0xFF33BAB70C0C34ABuL, 0x5BD7B3C8DDA2299FuL)),
                ((+1, 1, 0x87D7924DDD929390uL, 0xB55F61C3E0B7C99EuL), (+1, -4, 0xCCF4388B52C62DB0uL, 0xAE5A963E8B67CC8EuL)),
                ((+1, -1, 0x88F4615F748B568EuL, 0x902CCA97B038AFD2uL), (+1, -3, 0xDB3D54F3DAA6FB07uL, 0xC53BE33D4C6D09C9uL)),
                ((-1, -3, 0xF748AFF135FD3BA8uL, 0xCF51006CA743C351uL), (-1, -7, 0x980DC416C2CABB9FuL, 0xCB59B133EB53FF30uL)),
                ((-1, -7, 0xFD5F11146817BEAEuL, 0x290BD29DC6F194E8uL), 0),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p8125_0p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8A5CC00BBF73C050uL, 0x1D3B6BB43D92FD23uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0x8468B6C6A0EE758BuL, 0xC20FC43D7DAF1E61uL), (+1, 1, 0xB008923352EBF407uL, 0xD5224D260581E74FuL)),
                ((+1, 1, 0xBED5928618E0848FuL, 0xD6098B88C205F60AuL), (-1, -2, 0xC6523DEEF1773F33uL, 0xF94790D09275251CuL)),
                ((-1, 1, 0xD370EBADBCD60A12uL, 0xCD7F6B6AC7B00290uL), (-1, 2, 0xB00924C27B31114BuL, 0x840A1B47B7FEB89FuL)),
                ((-1, 2, 0xE5FCD3557C0C537AuL, 0x0D59110075A35628uL), (-1, 0, 0xD77914A9EA751815uL, 0xA62C19961D136D1DuL)),
                ((+1, 0, 0x82DE1E0A9321E1E8uL, 0x040A3803FB3F0CCCuL), (+1, 1, 0xC23C50A2EAAD8C7DuL, 0xD4D18EC7D7BF40ACuL)),
                ((+1, 2, 0x8E747AAD3CC936D5uL, 0x00331893D7B333DEuL), (+1, -1, 0xF7EE9D7B129F6076uL, 0x9EF6C75D6D1A99C3uL)),
                ((+1, -5, 0xD2938FA06A70B24FuL, 0x2DBA59C333338D47uL), (-1, -2, 0xCD405F04A1D1BBAFuL, 0xC8B30FCEACC5BDCFuL)),
                ((-1, -1, 0xA9795402F435AB03uL, 0xB204B0C0DD068DBEuL), (-1, -5, 0xE926AFDCDE0C076FuL, 0xE36F7F4ADC31711BuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p875_0p8125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xBD959266EDFA95FAuL, 0x1BBCB81A288ACA43uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xF93089A5016C3713uL, 0xBE40904C95320573uL), (+1, 1, 0xBD7DF4FEC8DA1107uL, 0xBE1FE129CC020ACFuL)),
                ((+1, 1, 0xE41F62B6A4890E43uL, 0x36419BE5BC94475DuL), (+1, -3, 0xC1DCC96F0DE0F5B1uL, 0x2C781E30960A4ED8uL)),
                ((-1, 1, 0x8C0D6AEBDF5C2726uL, 0xD3170961FD45F8F8uL), (-1, 2, 0xA71D9A3FEF1F1F81uL, 0xEBFEE2F60CF48A74uL)),
                ((-1, 2, 0xE88AE76B3FDFE2BDuL, 0xA7E3A79A3BE3FE06uL), (-1, 1, 0x8667F3C135049EDFuL, 0xED2ECAF3D54ADC94uL)),
                ((+1, -6, 0xFE6CB6F7CE021E99uL, 0xA3BF2AF9737F626CuL), (+1, 1, 0xA6F203697E352AA1uL, 0x4CEF5471A7F137AFuL)),
                ((+1, 2, 0x822BECC9AC0DCF04uL, 0xDCA540FDFE6DF846uL), (+1, -1, 0xEF7EB3D15449454CuL, 0xE32BE0DB852394B0uL)),
                ((+1, -3, 0xD013C38BAA7DEE18uL, 0x1216B30788920BEAuL), (-1, -2, 0xA6BE3B0E07525BE6uL, 0xD90CB8EC70D444BDuL)),
                ((-1, -1, 0x8EEE27B414E04BFFuL, 0xBB3DF1D68ADE7A47uL), (-1, -5, 0xBB90CC7650D58A23uL, 0xECC9D7B975235B60uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_1_0p875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                (0, (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xC90FDAA22168C234uL, 0xC4C6628B80D85100uL), (+1, 2, 0xD3548409D8EDED0AuL, 0xB0C9BA5C61CAABA4uL)),
                ((+1, 3, 0x996974712550498DuL, 0x632F02B59FECCC99uL), (+1, 3, 0xCC6680EA22C74C5FuL, 0xC55ACBD87D1D7071uL)),
                ((+1, 3, 0xF86AECF9AFDC6A04uL, 0x8204D95A06495D9BuL), (-1, 0, 0xD40191FFCAA43F3EuL, 0x6022EFD7698A7CA7uL)),
                ((-1, 3, 0x956909D4E6199D2EuL, 0xFF092E2E3EF51A2EuL), (-1, 4, 0xDADBD4F48F178F0BuL, 0x0D7BC25E5C95A8F8uL)),
                ((-1, 5, 0x993169F316A1EB7AuL, 0xF86110D52F36805AuL), (-1, 3, 0xD5C5EB552A543039uL, 0x8A03969F84E265C0uL)),
                ((-1, 2, 0xB281F954035040FFuL, 0x0C8D6DED353E98EDuL), (+1, 4, 0x97C4758291146E2AuL, 0xD3C60D4263E4F72DuL)),
                ((+1, 4, 0xF7488FB1294E7518uL, 0x185D7C32C405683CuL), (+1, 3, 0xB12E88FC7736781BuL, 0x785FF59F534B3BAAuL)),
                ((+1, 2, 0xD6B01E31E74893A0uL, 0xC9EC317353DFA836uL), (-1, 2, 0xAB4BC13ADF3164E8uL, 0x492B6EE901AC5BC9uL)),
                ((-1, 3, 0x9F855597377FEC18uL, 0x6CCA5971ED8ADA16uL), (-1, 1, 0x923538E03E7B448CuL, 0xE7FC795249C6017EuL)),
                ((-1, 0, 0x894FACA5643BD5ACuL, 0xD98877FED5F03F72uL), (+1, -2, 0xF5221AC95AD2FE1EuL, 0xA17224884C5314B8uL)),
                ((+1, -1, 0xE79091E4AE108679uL, 0xF47C36C7C26A0644uL), (+1, -5, 0xF60F72732096774CuL, 0x8B9FAEA95C837165uL)),
            }));

            public static ddouble Value(ddouble x) {
                Debug.Assert(x >= 0 && x <= 1);

                ddouble y;
                if (x <= 9.765625e-4) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_0_0p0009765625);
                }
                else if (x <= 0.001953125) {
                    y = ApproxUtil.Pade(x - 9.765625e-4, pade_0p0009765625_0p001953125);
                }
                else if (x <= 0.00390625) {
                    y = ApproxUtil.Pade(x - 0.001953125, pade_0p001953125_0p00390625);
                }
                else if (x <= 0.0078125) {
                    y = ApproxUtil.Pade(x - 0.00390625, pade_0p00390625_0p0078125);
                }
                else if (x <= 0.015625) {
                    y = ApproxUtil.Pade(x - 0.0078125, pade_0p0078125_0p015625);
                }
                else if (x <= 0.03125) {
                    y = ApproxUtil.Pade(x - 0.015625, pade_0p015625_0p03125);
                }
                else if (x <= 0.0625) {
                    y = ApproxUtil.Pade(x - 0.03125, pade_0p03125_0p0625);
                }
                else if (x <= 0.125) {
                    y = ApproxUtil.Pade(x - 0.0625, pade_0p0625_0p125);
                }
                else if (x <= 0.25) {
                    y = ApproxUtil.Pade(x - 0.125, pade_0p125_0p25);
                }
                else if (x <= 0.375) {
                    y = ApproxUtil.Pade(x - 0.25, pade_0p25_0p375);
                }
                else if (x <= 0.5) {
                    y = ApproxUtil.Pade(x - 0.375, pade_0p375_0p5);
                }
                else if (x <= 0.625) {
                    y = ApproxUtil.Pade(x - 0.5, pade_0p5_0p625);
                }
                else if (x <= 0.75) {
                    y = ApproxUtil.Pade(x - 0.625, pade_0p625_0p75);
                }
                else if (x <= 0.8125) {
                    y = ApproxUtil.Pade(0.8125 - x, pade_0p8125_0p75);
                }
                else if (x <= 0.875) {
                    y = ApproxUtil.Pade(0.875 - x, pade_0p875_0p8125);
                }
                else {
                    y = ApproxUtil.Pade(1d - x, pade_1_0p875);
                }

                return y;
            }
        };
    }
}
