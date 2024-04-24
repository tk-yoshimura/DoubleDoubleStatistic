using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class MapAiryDistribution : StableDistribution<MapAiryDistribution>,
        IAdditionOperators<MapAiryDistribution, MapAiryDistribution, MapAiryDistribution>,
        ISubtractionOperators<MapAiryDistribution, MapAiryDistribution, MapAiryDistribution>,
        IAdditionOperators<MapAiryDistribution, ddouble, MapAiryDistribution>,
        ISubtractionOperators<MapAiryDistribution, ddouble, MapAiryDistribution>,
        IMultiplyOperators<MapAiryDistribution, ddouble, MapAiryDistribution> {

        public override ddouble Mu { get; }

        public override ddouble C { get; }

        private readonly ddouble c_inv;

        private static readonly ddouble mode_base = "-1.1615872711359706852500000803029112987";
        private static readonly ddouble median_base = "-0.71671068545502205331700196278067230944440";
        private static readonly ddouble entropy_base = "2.0072768184106563460003025875575283708";

        public MapAiryDistribution() : this(mu: 0, c: 1) { }

        public MapAiryDistribution(ddouble mu, ddouble c) {
            ValidateLocation(mu);
            ValidateScale(c);

            Mu = mu;
            C = c;

            c_inv = 1d / c;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * c_inv;

            ddouble pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * c_inv;

            ddouble cdf = CDFPade.Value(u, interval != Interval.Lower);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble x = Mu + C * QuantilePade.Value(p, interval != Interval.Lower);

            return x;
        }

        public override ddouble Median => median_base * C + Mu;

        public override ddouble Mode => mode_base * C + Mu;

        public override ddouble Mean => Mu;
                
        public override ddouble Variance => NaN;
        
        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => entropy_base + Log(C);

        public override ddouble Alpha => 1.5d;

        public override ddouble Beta => 1d;

        public static MapAiryDistribution operator +(MapAiryDistribution dist1, MapAiryDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static MapAiryDistribution operator -(MapAiryDistribution dist1, MapAiryDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static MapAiryDistribution operator +(MapAiryDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.C);
        }

        public static MapAiryDistribution operator -(MapAiryDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.C);
        }

        public static MapAiryDistribution operator *(MapAiryDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public override string ToString() {
            return $"{typeof(MapAiryDistribution).Name}[mu={Mu},c={C}]";
        }

        private static class PDFPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xCA41ADEA0F2905FCuL, 0x536FE45EF780E1ADuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0xF4D5653B2AF6881CuL, 0xF388720EA9781DB1uL), (+1, 0, 0x90781480ED0964B9uL, 0x66BA8184A616E8B6uL)),
                ((+1, -5, 0xA3B7F98B417A2505uL, 0xE09D2344D00EADDAuL), (+1, -1, 0xCB09D758EF57D08BuL, 0xF0E4EA636758DF70uL)),
                ((+1, -7, 0xC92E7037324C277FuL, 0x713CA04384C62751uL), (+1, -2, 0xC4F78E98110B76B7uL, 0x528035223E944BE3uL)),
                ((+1, -9, 0xBE7843FBDC0F3C5DuL, 0x9E47C739B108CDC7uL), (+1, -3, 0x94FB2FC3E9300744uL, 0x50BA757297B9BCFFuL)),
                ((+1, -11, 0x85B09C4941C5F705uL, 0xA9800EDFB30988D3uL), (+1, -5, 0xB2866EAADBF6EC28uL, 0xAA81587CDE8D09BCuL)),
                ((+1, -14, 0xA8E6C97C58804E5EuL, 0x797F1EFAC8164137uL), (+1, -7, 0xAFDEA0D6C3132937uL, 0x7885837004526EA5uL)),
                ((+1, -17, 0xA1599F6ACDE0FF79uL, 0x03A9C411B1AC099AuL), (+1, -9, 0x8D481CBAEC253916uL, 0x5D63E3A92649FB51uL)),
                ((+1, -21, 0xD069FB1DE9C21ABDuL, 0x02EE35A1ED12B8B0uL), (+1, -12, 0xBB4FD1CDE8B7B05BuL, 0x1026A393CBE439C7uL)),
                ((+1, -25, 0xF9CC63639D09CB22uL, 0xB97582B4E8054FE9uL), (+1, -15, 0xC5E1491A69A2881BuL, 0x3DC4D155098758A9uL)),
                ((+1, -31, 0xD1BE32CCD7E0C3BFuL, 0x3DC2D73B62150609uL), (+1, -18, 0xA537E57629159915uL, 0x97216B01EB810AB1uL)),
                (Zero, (+1, -22, 0xC0C20D618B8BEF62uL, 0x507EE49035FDDA2BuL)),
                (Zero, (+1, -26, 0x9278B89D1D3E02E9uL, 0xF6BE13B07C9C6AD0uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xD99A406F35AD0588uL, 0x4B9995763BCB0B32uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xC3F1C5DB3F4A0DE7uL, 0xF067D936DE3F5683uL), (+1, 0, 0x90C5E4A8D64614BCuL, 0xD7FAD2B7D64E6DCBuL)),
                ((+1, -6, 0x82EB7445A3BE6EBCuL, 0xA321800A72C1AA12uL), (+1, -1, 0xB7E234596DC91385uL, 0x869F467275C7FD93uL)),
                ((+1, -8, 0x8BC19CF036EF9ED6uL, 0x0598D028FE9BBB7DuL), (+1, -2, 0x9CE5E7E9CB4E7DE7uL, 0xC450A7ECAC274A8DuL)),
                ((+1, -11, 0xBBE2DAE168607528uL, 0x5C5B80C07D30E010uL), (+1, -4, 0xC69D7B1E4499759BuL, 0x00F6A1EE73E9E243uL)),
                ((+1, -14, 0xD6E010EDFDAF2F9FuL, 0x54CA2CBA4EA64A25uL), (+1, -6, 0xBF9FE9D5F8031983uL, 0xC1E4523705095663uL)),
                ((+1, -17, 0xAE91CB93089F94CCuL, 0x54DF3B8E4DDF4826uL), (+1, -8, 0x8F81258C8F24B692uL, 0x4277ED855263AD02uL)),
                ((+1, -22, 0xF139CA86015B6BD7uL, 0xFD0939E7ED4C3232uL), (+1, -11, 0xA429796B43B167C6uL, 0x3383CF01AF8DBB9EuL)),
                ((+1, -25, 0x998B4A93DB5DC046uL, 0xDF6B2E3A32A601C8uL), (+1, -14, 0x8BEFA5D9CD98BBD5uL, 0x7AD83AE256606530uL)),
                ((-1, -29, 0x9180F79B09802077uL, 0x78B687A61839D7E1uL), (+1, -18, 0xA1400A63347CA5C1uL, 0xD8F53799495C5151uL)),
                ((+1, -34, 0xE376F49061D4CB4DuL, 0xF177E69094D4B573uL), (+1, -23, 0xD299BCCAC2FFB5ABuL, 0x2936BB47302908D5uL)),
                ((-1, -39, 0xAE633639CB0F69D7uL, 0x2190989C9841ADC1uL), (-1, -33, 0x873D362B9379F61DuL, 0x38B682AFC6C797E5uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xDAA971E447E1363FuL, 0x3CB7A4FB4CB38B64uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xA728F5855A4B289DuL, 0x66B1341028573870uL), (+1, 0, 0xB89C4F18280B9953uL, 0xFB060E1B06E4E2ACuL)),
                ((+1, -6, 0x9DD5B6EA0B2B29A8uL, 0x9F92395911494C93uL), (+1, 0, 0x8AC63BA430888489uL, 0x5AD93E5B2C987B33uL)),
                ((+1, -8, 0xCCC4180A324819F5uL, 0xDAB4C6F7CD7BAB3EuL), (+1, -1, 0x8B09D54EE261E119uL, 0xD82D77B1B821FB1FuL)),
                ((+1, -10, 0xC64C8D140960548CuL, 0x1FF0D280F2FE5423uL), (+1, -3, 0xCCC88E00166B9071uL, 0x77335B7BD0541E24uL)),
                ((+1, -12, 0x92E361054F0ACBDBuL, 0x9C423577E5881771uL), (+1, -5, 0xE80B711377D6DF4AuL, 0x66378EAE4FD594E8uL)),
                ((+1, -15, 0xA5EE006AE81CD839uL, 0x0D44A3F784A047DBuL), (+1, -7, 0xCEF7F72040FBBBE9uL, 0x0295E115999898E1uL)),
                ((+1, -18, 0x8CFC8EC62BA3CE17uL, 0x66DBE52F552AE4C6uL), (+1, -9, 0x924FD40CC2E0A5D4uL, 0x7F6CBA54ACC92AB6uL)),
                ((+1, -22, 0xA7BBD9883F053902uL, 0x82F3E41ABFADB411uL), (+1, -12, 0xA30D10D043577F3FuL, 0xAB42FF467FEF9E09uL)),
                ((+1, -27, 0xFC31C21359B03BA9uL, 0x04CF9436BB3C699DuL), (+1, -15, 0x8C2C2CA69445F0AFuL, 0xC6967F361D48768AuL)),
                ((+1, -33, 0xFC22111B0A520346uL, 0x679F9279750FC80FuL), (+1, -19, 0xB1DE265D3051E4E5uL, 0xFD1F1E6A64ADF019uL)),
                ((-1, -39, 0xAF5E94CA38682FC8uL, 0xA9CE691D47B0D0F0uL), (+1, -23, 0x973CDE9310D1CD42uL, 0x44DCE616D2AD922CuL)),
                ((+1, -45, 0xC01E5CEE889CF09DuL, 0xF3B2076EEBEA9400uL), (+1, -28, 0x86DC25C24D9B84D0uL, 0x5E1A8601F88A9291uL)),
                ((-1, -51, 0x8819E2F3CA63670CuL, 0xFDF72EC8F7178428uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0x8236580FAAEBE5ADuL, 0xAD68B145ED634949uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -7, 0xFB313EA20A20B93AuL, 0x508AA27B7A92AC4AuL), (+1, 0, 0xBEAA49812E2E3D78uL, 0x03E8749168194979uL)),
                ((+1, -8, 0xF9D4DDA07A688039uL, 0x9279F68DC4E7E2F1uL), (+1, 0, 0x8ABFEE7B96C20AF1uL, 0x24D3FE618CCD5FE8uL)),
                ((+1, -9, 0xA1C1C451D0CB3F65uL, 0x6F71934146E87CB1uL), (+1, -1, 0x81A4470FB28842A1uL, 0x0620C2E6E39D8151uL)),
                ((+1, -11, 0x9516D8B1D82591D7uL, 0xA034FB2DB72913E9uL), (+1, -3, 0xACBBC9E29F331AD4uL, 0xEAEDAEF84CFC370AuL)),
                ((+1, -14, 0xCAEBE390054D28F4uL, 0x6DBCA23543382203uL), (+1, -5, 0xACB94576AF19EEE4uL, 0x0667A37D3EFF51F5uL)),
                ((+1, -17, 0xCE0C5709833FFBF3uL, 0xE0BCBAADBC2B4160uL), (+1, -7, 0x8506BD1AE9EBB52DuL, 0xE10C964042D0D940uL)),
                ((+1, -20, 0x9A8A48229A59A9C7uL, 0x1C9C3992B16A5F90uL), (+1, -10, 0x9F7A391199B1EB63uL, 0x7C9C4C173ABFCD5EuL)),
                ((+1, -24, 0xA59C3A71F9D3384AuL, 0x2E60DB79F3A172B5uL), (+1, -13, 0x948EE61883F7A65EuL, 0xBD2897F7DEF6A928uL)),
                ((+1, -29, 0xED14AFB8B6E02416uL, 0xC280101F8B48864DuL), (+1, -17, 0xD4225D0269512D9DuL, 0x150115A80E755802uL)),
                ((+1, -34, 0xC53E5B3C53E6865AuL, 0xC5AE0107268E3FE6uL), (+1, -21, 0xE19275F25C44B0D2uL, 0x68216358B31A64BCuL)),
                ((+1, -40, 0x8E85AB0B83D8EDFBuL, 0x4DDA4D1F2256C3B9uL), (+1, -25, 0xA98F9713EEA4A7FBuL, 0x3F18608C2F551DECuL)),
                ((+1, -49, 0xA62A933C9DE11D54uL, 0x048FF180A7021518uL), (+1, -30, 0xA401F62362DBC4E0uL, 0x56FBAA6BBB323BEBuL)),
                ((-1, -59, 0xA960E4E7960E3896uL, 0x64F2709B6A240B58uL), (+1, -36, 0xA97BBDCFAF55A222uL, 0x1B24F91014FE46CBuL)),
                (Zero, (+1, -44, 0xF628E2630BA759F0uL, 0xE42E194B38149836uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -9, 0xD35D775CA2E3F28DuL, 0x44D071481C968757uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -9, 0x9E1170C1E3552A20uL, 0xA2566567D213ED6BuL), (+1, 0, 0x8694C74C70AB2BD7uL, 0x0195C4C0B25E1349uL)),
                ((+1, -11, 0xD8B474A299C7527BuL, 0xBB91E70FC3569C9BuL), (+1, -1, 0x83129D1A9A2D3A80uL, 0xBE391BC6D5BD2143uL)),
                ((+1, -13, 0xB100306324DF6B7BuL, 0x101BF7C81F1ED1F9uL), (+1, -3, 0x9BD47D642B4A9B02uL, 0x6C570ABA80CB8EDEuL)),
                ((+1, -16, 0xBC928676337720FEuL, 0xF722D45342DB8766uL), (+1, -6, 0xFAFF32081A47367CuL, 0x9E9F2086E4074D49uL)),
                ((+1, -19, 0x871FA0E533ECFDEEuL, 0x9CEFB905B696F3E7uL), (+1, -8, 0x8FB5249AC14A858DuL, 0x297E116338ACD0F1uL)),
                ((+1, -23, 0x81D1117117B75534uL, 0x43BD6FF666AE83DCuL), (+1, -12, 0xEEE3E96C1C9C965AuL, 0x5BBEEE076182FFF3uL)),
                ((+1, -28, 0xA1D6B925A87D6569uL, 0xF64972878958DED0uL), (+1, -15, 0x9082DE87B3C2C65FuL, 0xA323F26E64EE0D90uL)),
                ((+1, -34, 0xF50554CAD906F9C4uL, 0xEF5C92AF753F5E3BuL), (+1, -20, 0xFB40DEBFEBCE970FuL, 0xA13ADEB17F2027E6uL)),
                ((+1, -40, 0xC8CB010A25730E0AuL, 0x3C23F8B24CF9A907uL), (+1, -24, 0x985854B534721F6EuL, 0x6E243A9FFD7A5205uL)),
                ((+1, -47, 0x8F2F859253095485uL, 0xEA03422AA2B1FFA2uL), (+1, -30, 0xF49CF4551841C23DuL, 0xED150CBCADE1597BuL)),
                ((+1, -57, 0xB8C18EDCEA889F2EuL, 0x2FA034A65CE1589AuL), (+1, -36, 0xEE45E876AFE50A8DuL, 0xD6C6865E036265E2uL)),
                ((-1, -67, 0x8FD3E28E09283D67uL, 0x212A0668A459CC38uL), (+1, -43, 0xF135C70E8CC9E87FuL, 0x5E0BAD7DC838DEAFuL)),
                ((+1, -77, 0x984BF6359C8C620EuL, 0x7F8A1EEF655DC7F5uL), (+1, -51, 0xB3C7AA58B82F32C8uL, 0x13B93A327D409D4EuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -11, 0x98B4C50CE7201D4FuL, 0x325FFE0BF8352867uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -13, 0xE8ADA879B9145346uL, 0xF1854245B99E833DuL), (+1, -1, 0x895D4FCCD322781EuL, 0xD10AA725F14FC717uL)),
                ((+1, -15, 0x9A051BBA7B3EEA60uL, 0x149724CB8D44E568uL), (+1, -3, 0x84BF453C70FBB960uL, 0xC957E5F5CFC05AD4uL)),
                ((+1, -19, 0xE7665DB438E1C7FBuL, 0xDFF31EC5E87F0762uL), (+1, -6, 0x984A08D63DC7CE03uL, 0x29191E37B427F40DuL)),
                ((+1, -23, 0xD86D8A13E876DE95uL, 0x22134F0F70DAC3EAuL), (+1, -10, 0xE62FFCA4849EF871uL, 0xA78CEAB4629CDA73uL)),
                ((+1, -27, 0x824B25402FE23843uL, 0x39EEF0627928C0CFuL), (+1, -14, 0xF0855136E558CC49uL, 0x79F0FDC8E740962FuL)),
                ((+1, -33, 0xCA2E82B2B99A3DCAuL, 0xE39AF7765FD2C7BDuL), (+1, -18, 0xB15CE03CC16C031DuL, 0x12AFF263D5D0AFAEuL)),
                ((+1, -39, 0xC547488600AEE971uL, 0xB83A6D4DACF7564FuL), (+1, -23, 0xB92271B8E9D2943DuL, 0x9605532CE5CE7EE2uL)),
                ((+1, -46, 0xE5686ACB0FA3CA46uL, 0xCEBEDC2D87F75FFDuL), (+1, -28, 0x87248B1A842959F7uL, 0x4A6348348219467AuL)),
                ((+1, -53, 0x8FC1D19F4CDA0928uL, 0x4D63805EE9E19E00uL), (+1, -34, 0x86355AA5FF1DF115uL, 0xAB244066FEAB5CABuL)),
                ((+1, -62, 0x9E18591A0DAF2DCEuL, 0x965F4B3CAE5094BBuL), (+1, -41, 0xACC755719AD79264uL, 0x7F2D3048C93C7215uL)),
                ((+1, -73, 0xA015190BFBF4B6D2uL, 0xFC31E8751ED4CA6FuL), (+1, -48, 0x84F4C6A9B468D34FuL, 0x87236076865CA98FuL)),
                ((-1, -85, 0xC81AECA702F88436uL, 0xE1FCF3C4BB4A2282uL), (+1, -57, 0xD350859F55BB75CFuL, 0x6276F812C813E14FuL)),
                ((+1, -96, 0xAE6C886A69AD5E1CuL, 0xB5334366CB371F17uL), (+1, -67, 0xF807E2B006584F85uL, 0x4848D59CF2D89CFEuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -14, 0xD88FE560FF6B053DuL, 0xC0D11DFC08376236uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -16, 0x824BA414282A0581uL, 0x462D201ECBBEE51AuL), (+1, -3, 0xE9FC2BF5F18F9700uL, 0xCCBF57B6113955ACuL)),
                ((+1, -20, 0x844C34157564B084uL, 0x001651F5D466D898uL), (+1, -6, 0xBD68A4C67413E25DuL, 0xBAEF30377EF10E26uL)),
                ((+1, -25, 0x93B48BB20905A6A7uL, 0xFE5CDD5E748F8D24uL), (+1, -10, 0xB2AF419912F73A42uL, 0x3007A4543B2F1CE5uL)),
                ((+1, -31, 0xC6442DAB0F8B6789uL, 0x36C3E69543617472uL), (+1, -15, 0xD9830C4551F85244uL, 0x8136ADC9C93D93EEuL)),
                ((+1, -37, 0xA47CD8242CB9E73CuL, 0xA1CF0C8F04742177uL), (+1, -20, 0xB2AEE2C755252CEFuL, 0xE87F136E3257836DuL)),
                ((+1, -44, 0xA72C44B569C0C8D9uL, 0x90884A4DC79A20AFuL), (+1, -26, 0xC95C55A72FEBCA7FuL, 0xA9BF531EAA536829uL)),
                ((+1, -52, 0xC73B6A6FFD6DA26EuL, 0xECCF337B7C524082uL), (+1, -32, 0x9B18B35BC8E87010uL, 0x37C9936DEA7CABADuL)),
                ((+1, -61, 0xFD8DADB4EA82D59AuL, 0x4910C343831F75F0uL), (+1, -39, 0x9FCD092915DAE036uL, 0x78D678018D2A8159uL)),
                ((+1, -70, 0x8CF9079EC9A5162DuL, 0xD8A5FF98E21B32EFuL), (+1, -47, 0xD2DB6F4024DCDBB4uL, 0xD1ABE8E185A59E3EuL)),
                ((+1, -82, 0x903AA0C9F228BB7CuL, 0xA7AD43B8FB997F39uL), (+1, -55, 0xA4F0B5F3F9FF150FuL, 0x97232D1B6E972758uL)),
                ((-1, -95, 0xB67BBB3B80965A7AuL, 0xECEEBEAC929A1065uL), (+1, -64, 0x849AB8E7F9A637E9uL, 0xA3916248FF0DD1B5uL)),
                ((+1, -107, 0xA18D1033F91DEB5FuL, 0x3CFC1D55B85130B6uL), (+1, -75, 0x9D2FD170846F3AC2uL, 0x9617A9B18E7F3BFCuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x99319F36F945E34EuL, 0x5EC369E04605551FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -11, 0xE262A7B9FA58263EuL, 0x28F4D1E8DD0C5A13uL), (+1, -10, 0xBD2797AA3F557413uL, 0xFBC72B25B6FA7B20uL)),
                ((+1, 6, 0xB15F468AAA704A76uL, 0xAE26E768E8C580F6uL), (+1, 7, 0xA153B862DF774E6CuL, 0xBF2FBB2F4F1C642FuL)),
                ((+1, -5, 0xEEE91E0010D0F7DCuL, 0x80851C41E7BFB891uL), (+1, -4, 0xEE69591CB9EAA866uL, 0xCE2D61CBDAFA3355uL)),
                ((+1, 11, 0xB0E10D2818E95B21uL, 0xA465E445EDB7EDE6uL), (+1, 12, 0xCAF647CC6C362E4BuL, 0x743AD3CBB7E69924uL)),
                ((+1, -2, 0xA522905A0DB0AA79uL, 0x1AFAC2A40C3C6AECuL), (+1, 0, 0xC78DF9A53A372B30uL, 0xDA278C939D6AD9DAuL)),
                ((+1, 13, 0x8DDF6F2C221DBBE8uL, 0x91532539B0778434uL), (+1, 15, 0xE1FC8D9082955721uL, 0x8306E0B6F127BFD8uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_0p5_0 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xFC6E98982981FE2CuL, 0xC0A373B5A7059DCAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x87FDD730E8602927uL, 0x6ACD7A56AF602363uL), (+1, -1, 0xE3B85B0D6FF847FDuL, 0xA32AC301ECB4525CuL)),
                ((+1, -6, 0xE8A67BC9371B5473uL, 0x6694BCDE5FAEFFEDuL), (+1, -1, 0x936D9CCF77935B60uL, 0xC4B64B6D0566292CuL)),
                ((+1, -8, 0xF003E651C0A3BD10uL, 0x8DFFF87EE04795D1uL), (+1, -3, 0xFAC1AA74FCB27123uL, 0x5D25B0A6BC766B0EuL)),
                ((+1, -10, 0xF10F33FC71ADFFB9uL, 0x229CF9CE6832616AuL), (+1, -4, 0xAE73206872796DF5uL, 0xB354AA630CB603ACuL)),
                ((+1, -13, 0xC6909D0582632659uL, 0x84E9A5B97BAD0B13uL), (+1, -6, 0xB7FECAFB23E446E7uL, 0x819EE14B3CED7EF7uL)),
                ((+1, -16, 0xAF7FC21E0FA090EAuL, 0x340A5D067422BEA1uL), (+1, -8, 0xA36F5303B1F4B3C1uL, 0xB71A3908399488D5uL)),
                ((+1, -19, 0xECA3F6032F56E098uL, 0xD9A13F20EB444670uL), (+1, -11, 0xDC803F5F6BFEC470uL, 0x2C0347FF0353C6E3uL)),
                ((-1, -24, 0xD86EB4C2B0E5F419uL, 0x1EA2ADA702B8155FuL), (+1, -14, 0xF86ECC3585A812E7uL, 0xDE1376FCFF3272DEuL)),
                ((+1, -28, 0xD362B796255008ABuL, 0xB50A40975AF40D35uL), (+1, -17, 0xB81C1E278811326CuL, 0x54CA5C04B58E1488uL)),
                (Zero, (+1, -21, 0xD73B732F6A1E81BAuL, 0xBDCE1EDE36CC140AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8DC093A4A1102639uL, 0x1975F0183CC69FFCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xBC45DB5FEE8C3F8BuL, 0x3CDF1E3760A309BDuL), (+1, -1, 0xC380B8DAF0F8F146uL, 0xB2E21D4341EB9931uL)),
                ((+1, -5, 0x97825FEEC75F3806uL, 0xA2A8F5C333CC7F0BuL), (+1, -2, 0xFF3BE52148BD1D73uL, 0x7C1FF04AD57A9218uL)),
                ((+1, -8, 0xA391249D3352E6E2uL, 0x58B1CBB8D9923E24uL), (+1, -3, 0xCD46EAD5FCF331D2uL, 0xF061266C30E34F7AuL)),
                ((+1, -9, 0x88074EF687C4B3E1uL, 0x405087E2AB8258AAuL), (+1, -4, 0x91D80681D7CEDD27uL, 0xF4F0B948353A8A84uL)),
                ((+1, -12, 0xAE4DD3F1A89A89C2uL, 0xE6FEC959A580DAFBuL), (+1, -6, 0x964115B66FB2A07CuL, 0xAF1F18ED1C5C03FAuL)),
                ((-1, -18, 0xF1BADB320DF1C670uL, 0x3AF93BE426FDB474uL), (+1, -8, 0x8B4EE130A5A5F6CEuL, 0xD7B99410BE089ACAuL)),
                ((+1, -19, 0xE1D562EE19EACCBDuL, 0x8A65BC4197EC9212uL), (+1, -11, 0xBABAF018B10BFDFCuL, 0x3A57ADFDA040128AuL)),
                ((+1, -21, 0xC9DFE45C08598DDEuL, 0x56BE6004C5F2BBEBuL), (+1, -14, 0xE469E737B16D489BuL, 0xB7388BA18182FB89uL)),
                ((-1, -24, 0xF5B83C472E999D5FuL, 0x38C55BE7E8BD49CCuL), (+1, -17, 0xA8EE2AA38796B3EAuL, 0x967CA0DD447EEB03uL)),
                ((+1, -28, 0xCCB1EB9E55668C34uL, 0xE6A1CB5BADBC9F06uL), (+1, -21, 0xF0273259489C876CuL, 0x54B5CDE25441B526uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xDBA1A65E27582E2AuL, 0x573BE75E924E6455uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x847539359BFF4643uL, 0x1311366213457FFEuL), (+1, -1, 0x886B90F56DB02C56uL, 0x9120F66E069B37ACuL)),
                ((+1, -4, 0xEA339CC0B6922F68uL, 0x1D7A3632640243F2uL), (+1, -2, 0xDC83C069585D7834uL, 0x4B3567F561F922F3uL)),
                ((+1, -6, 0xAA85193DF3544844uL, 0x1AB6542E59FBD9B7uL), (+1, -3, 0x9B128C0266C7A9E2uL, 0x74E09B678AA35C1CuL)),
                ((+1, -10, 0xE9A9B90CDE546832uL, 0x97BD3587E56A50FFuL), (+1, -4, 0x8BDFAA6A5EE3408FuL, 0x00F5C98C0931A9D2uL)),
                ((+1, -11, 0xBAD00AEFBEF9CCB7uL, 0x8C230AC21EFC0F32uL), (+1, -6, 0x8C96D0B311B64D2FuL, 0x7DD94C8BB1CD9E2FuL)),
                ((+1, -13, 0xE38DAA7B57355B01uL, 0x2FDC64588063DD78uL), (+1, -8, 0xB28AE8AEDCB7CDA3uL, 0x497F0F0FF93F08D0uL)),
                ((+1, -17, 0x8AD6A89F2B4B8701uL, 0x618EDA76A90320DEuL), (+1, -10, 0x8111586902A7927DuL, 0xB8E360BD6738C25BuL)),
                ((-1, -20, 0xDA6540C0467C12C8uL, 0x6B78E4012E7DFEDFuL), (+1, -13, 0xFA32D43B7864A8B7uL, 0xC21D6D85534C0DD7uL)),
                ((+1, -21, 0x8B5549EBDFF37E84uL, 0x50D57EE879968395uL), (+1, -16, 0xF37B0E1C78D3303FuL, 0x780E3FEFBD2C7C1BuL)),
                ((+1, -25, 0xA345CC69D101AD09uL, 0x748E02DB80ED1668uL), (+1, -18, 0xBB26A911C00635BAuL, 0x95ABF95A4D462171uL)),
                ((-1, -27, 0xA4CBB0F395CDDC44uL, 0x911A4BD0DC8F393AuL), (+1, -22, 0xBF91E8137C046901uL, 0xD9D6A84C105729A7uL)),
                ((+1, -31, 0x942B271AACF16A42uL, 0x05FFB93AD08D3CCCuL), (+1, -25, 0xF03A5FFC9F1D19D5uL, 0xA0349DF26F223EC3uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8C7229BC422E08F1uL, 0xA02EBF480F68BE60uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xC3F2DD3DBBE1D77CuL, 0xD599A7608FA8C5BEuL), (+1, 1, 0xB4C939C8B800AC41uL, 0x56347518CA2FC6EFuL)),
                ((+1, -1, 0xFD18D0201B7259F4uL, 0xFC73DA285997EE60uL), (+1, 1, 0xEB50DE505D8FEE19uL, 0xAD3367639F70B9E1uL)),
                ((+1, -1, 0xCA73DDCE123913C0uL, 0x7641CE1A467440A9uL), (+1, 1, 0xBD122C8711B5E044uL, 0x83FDCD1E036299FAuL)),
                ((+1, -2, 0xE1E9210D8909692BuL, 0x44CCF8A227D05851uL), (+1, 0, 0xD385D0F63EBB955CuL, 0x60DB6EE47AA2D1BFuL)),
                ((+1, -3, 0xBB279F6A52498469uL, 0xAD2D35AE67B85A50uL), (+1, -1, 0xAF8253529B5C19D3uL, 0x46E5856C8F2006C9uL)),
                ((+1, -5, 0xEE6853AADB10F648uL, 0x0C19D2BCDFF9265BuL), (+1, -3, 0xDFC9E5CF65A79231uL, 0x7FEE4AF6C299BB7FuL)),
                ((+1, -7, 0xEDC0D69DF1654B46uL, 0xB184C3B81E0CD599uL), (+1, -5, 0xDF4F55E8082FDB6FuL, 0x5E671DD5887B2007uL)),
                ((+1, -9, 0xBA9DDB1FB414A1E5uL, 0xB799FBFA4B04D380uL), (+1, -7, 0xAF589A49AB23F37AuL, 0x138792543DB02616uL)),
                ((+1, -12, 0xE5008247553AE5B3uL, 0x40EA47A907CF767AuL), (+1, -10, 0xD7387B2D64674590uL, 0xBBB78625F0AEBFFCuL)),
                ((+1, -15, 0xD6DC08EB1E5234A2uL, 0x53CE5B8C4115E24BuL), (+1, -13, 0xC9F44EFA593E34C9uL, 0x03DF5C767A6B5AF0uL)),
                ((+1, -18, 0x9318BE9B2C79F0CCuL, 0x6C2EBE1132E63BA0uL), (+1, -16, 0x8A44CBE26F198846uL, 0x6634B416A42DC092uL)),
                ((+1, -22, 0x851D0F33C0FD7053uL, 0x0B3AED3CE7260087uL), (+1, -21, 0xFA3FE4584D5B9778uL, 0xB602652BAF607D59uL)),
                ((+1, -28, 0xF4AB6C886982945CuL, 0xD654EEFB00860F73uL), (+1, -26, 0xE5FC65588F467058uL, 0xBCDF3B2C400CDD3BuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x88E78B8142BD3618uL, 0x68F4C6A729C28540uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x94E818A9BD397CB0uL, 0xFC9B4CBCC0033331uL), (+1, 1, 0x8B75945F841716C9uL, 0x639F50C55BBF71EBuL)),
                ((+1, -1, 0x99C2E7C05E8D105DuL, 0xD1971C1B4BC13378uL), (+1, 1, 0x902BEEA2A3D90016uL, 0x2128C13952D82EE1uL)),
                ((+1, -2, 0xC7B587F7B5E54C44uL, 0x7BDF507576ABE84FuL), (+1, 0, 0xBB671349F4366F1BuL, 0x928AE15C3D36D7B3uL)),
                ((+1, -3, 0xB634D7F428606B14uL, 0xA8D62E3338B1DA10uL), (+1, -1, 0xAB12E59742669F43uL, 0x398761CCC9485271uL)),
                ((+1, -5, 0xF6F3ADEBB74206A4uL, 0x07A160DF77A58A3EuL), (+1, -3, 0xE7F46FED36BF3003uL, 0xC2E2D2BE36CA1B18uL)),
                ((+1, -6, 0x8012042DB1ADB84EuL, 0x7DA0EFFFD6AAE595uL), (+1, -5, 0xF0A71DC65DA8B5F6uL, 0x35F8FCFB75FF92E2uL)),
                ((+1, -9, 0xCE2E0D7A4B0B83B0uL, 0xCBE751D68A8C49D6uL), (+1, -7, 0xC1C02F730E8CF994uL, 0xC9D139DC307034F5uL)),
                ((+1, -11, 0x8106BB663FA3499BuL, 0x7D5C3086A32ADFF0uL), (+1, -10, 0xF2875127D6227189uL, 0xE0CEDAB9FDFE7859uL)),
                ((+1, -15, 0xF87CFD48834F8B85uL, 0x22D870E03795A92CuL), (+1, -13, 0xE98F3B7A3E268425uL, 0x6A67CA1D962E3941uL)),
                ((+1, -18, 0xB35B3A5D9028521CuL, 0xCA5A12DB304AD94FuL), (+1, -16, 0xA896CCAE31F9BC33uL, 0x3DCC669EAD062F04uL)),
                ((+1, -22, 0xB84CD2CFE8822B79uL, 0xC5D0A0C0BBB9BF85uL), (+1, -20, 0xAD3D4B78D6F98ECAuL, 0xEB7438D1D4E15F6CuL)),
                ((+1, -27, 0xF27F873C7A142E24uL, 0x9CAA1215B12A1EF8uL), (+1, -25, 0xE3F1DCE12D394F3EuL, 0xD991C6E37786B4FBuL)),
                ((+1, -32, 0x9B01C7EBEB544A5AuL, 0x7B3C9FA414CB4C37uL), (+1, -30, 0x91B44D838046D79FuL, 0x7AF35F7853EAA11EuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8845530100984E7AuL, 0xB48C8C08A9BA7412uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xC9480B92B2F1AE37uL, 0x5200FEDCFB9FA9FBuL), (+1, 0, 0xBD19626F75133E7FuL, 0xD070641A0E0B5696uL)),
                ((+1, -2, 0x855726380D9A2254uL, 0x2B858FE7AC0AC6F9uL), (+1, -1, 0xFA945401E924557BuL, 0x14884D4E2172401FuL)),
                ((+1, -4, 0xD2E96CC240616841uL, 0xB988F857E6091E32uL), (+1, -2, 0xC633BB0AFFD78797uL, 0xF94CFFCED4F6E9E1uL)),
                ((+1, -6, 0xDE631C451A7998F1uL, 0x0DE5E051837C817AuL), (+1, -4, 0xD1016D1FAF7A62B9uL, 0xCD0E6BA4B71447E9uL)),
                ((+1, -8, 0xA4A8F88C2AB8886DuL, 0x0C68C87923C938B6uL), (+1, -6, 0x9AC33F275F2C005FuL, 0x0C83DB77883371F0uL)),
                ((+1, -11, 0xAF18DAD9396E1C88uL, 0xFFC6DADDBE21C7A3uL), (+1, -9, 0xA4948786704EBB2EuL, 0x02E2B5836C28E503uL)),
                ((+1, -14, 0x85E6A87ACA5ACCA6uL, 0xD33AF5FB4A7C3316uL), (+1, -13, 0xFBB979B3A4EC1270uL, 0xCA7A3D0D74E5CF1CuL)),
                ((+1, -18, 0x90785821E665FEDDuL, 0x90CE2ECF7C768E74uL), (+1, -16, 0x87CC8C15EA605E96uL, 0x0D6C6EA571E035FBuL)),
                ((+1, -23, 0xD19218B86BFE9A73uL, 0x022D2F6BC5DAE2E6uL), (+1, -21, 0xC4FE5180A25E45C2uL, 0xF05C5E4F6834ABCDuL)),
                ((+1, -28, 0xB7FBE6EF70E8065DuL, 0x106E78F2B0597747uL), (+1, -26, 0xACF13AD89BF99677uL, 0xA42549EFABB5C9EFuL)),
                ((+1, -34, 0x941C38A3AA71976BuL, 0x4A06B71A0AFBAD1EuL), (+1, -32, 0x8B38B3BC4FC465ACuL, 0xE92F0BCBB66EEDA1uL)),
                ((-1, -90, 0x8F9AD41E23441198uL, 0x5A3977A943FE602EuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x882F4B7A1409F5EBuL, 0x5A91D6F30B26D80FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x8927827821A82E76uL, 0xDBA4271421139069uL), (+1, -1, 0x80EA6B1882410041uL, 0x92472A90F1DFDA42uL)),
                ((+1, -6, 0xF264FCB5B7843B06uL, 0xC7CDA673864A986AuL), (+1, -4, 0xE3D6FB698BCD4F38uL, 0xFDD8BBF90BA1B5DEuL)),
                ((+1, -9, 0xF569CFDEE20C1A72uL, 0x0F0B8D8EB3E647C7uL), (+1, -7, 0xE6AE5E133ED0C4CCuL, 0xCE40959EEE362909uL)),
                ((+1, -12, 0x9BA744B8BA163F8CuL, 0xD60E8F3BC66D09BAuL), (+1, -10, 0x924F9A63D9282D86uL, 0xAF9D55075E313C7DuL)),
                ((+1, -17, 0xFD4289924515AF7BuL, 0x652BB1031D3C96BBuL), (+1, -15, 0xEE0F6D05D0B7FA5AuL, 0x5C629BCF44528AECuL)),
                ((+1, -21, 0x81032877020D7BCDuL, 0xA1F85D52C840F596uL), (+1, -20, 0xF28A1B2E8CEA41B5uL, 0x218EA114643DA7F4uL)),
                ((+1, -27, 0x9675BC64C78E9696uL, 0x2F0B364083D9815CuL), (+1, -25, 0x8D6E1E01F9CC5504uL, 0xF84596AD86EC783FuL)),
                ((+1, -34, 0x99BFD2CDB83F02F0uL, 0x3124A11F0CACA24FuL), (+1, -32, 0x9085AAD9736C6736uL, 0x905079FB24DA5F26uL)),
            }));

            private static ddouble PlusValue(ddouble x) {
                Debug.Assert(x >= 0);

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
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    ddouble u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u / x;
                }

                return y;
            }

            private static ddouble MinusValue(ddouble x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 2d) {
                    ddouble y;
                    if (x <= 0.5d) {
                        y = ApproxUtil.Pade(0.5d - x, pade_minus_0p5_0);
                    }
                    else if (x <= 1d) {
                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0p5);
                    }
                    else {
                        y = ApproxUtil.Pade(2d - x, pade_minus_2_1);
                    }

                    return y;
                }
                else if (x <= 32d) {
                    ddouble v;
                    if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else if (x <= 8d) {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }
                    else if (x <= 16d) {
                        v = ApproxUtil.Pade(x - 8d, pade_minus_8_16);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 16d, pade_minus_16_32);
                    }

                    ddouble y = v * Sqrt(x) * Exp(-(2d * x * x * x) / 27d);

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
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xAAAAAAAAAAAAAAAAuL, 0xAAAAAAAAAAAAAAAAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x8CEAA9B5C8264533uL, 0x9400D45E497FF8DFuL), (+1, 0, 0x80B0A0DBF0BDBC31uL, 0xF6CA4546FA0726D1uL)),
                ((+1, -5, 0xFBC3F04019274571uL, 0xBD3746E398E597C2uL), (+1, -1, 0xA00478381A2F10B2uL, 0x8B18D89FECD069F2uL)),
                ((+1, -6, 0x8BD3A81263DF56C0uL, 0xD4CEFCF06E72B528uL), (+1, -2, 0x88074BB1853BF00FuL, 0x102541DA31694FADuL)),
                ((+1, -8, 0x8AF1FC59C24E8DDFuL, 0x120D683E15CA248AuL), (+1, -4, 0xB20AE3CC91006DFFuL, 0xD62D7598711CE0EDuL)),
                ((+1, -11, 0xB7C879B3B9C8DD61uL, 0x343B4CF276033764uL), (+1, -6, 0xB5C04F53AED912AFuL, 0x99CC7DA706681471uL)),
                ((+1, -14, 0xFB34894D2BF3898FuL, 0x167666BDAC06DD6EuL), (+1, -8, 0x954CC0FB64647EABuL, 0x9E94E6FA8E711278uL)),
                ((+1, -17, 0xB5486E82B4C103AAuL, 0x6C1932DAD9ECE4EEuL), (+1, -11, 0xC28FC5619CDD51DDuL, 0xC2F4C7A9A93885FCuL)),
                ((+1, -20, 0x9D00A7F31E6DA46FuL, 0xA266C4D537EB3F4BuL), (+1, -14, 0xC89C0F03E22CC96AuL, 0x30273B0E26D9D08BuL)),
                ((+1, -25, 0xB4ACA328EA3F33CEuL, 0xA7052A5FD351277FuL), (+1, -17, 0x9AEE8BBF651DCBA0uL, 0x70BA9F43C69C6196uL)),
                ((+1, -30, 0xC3CA1F92AF3448DFuL, 0x45B14264666D1A18uL), (+1, -21, 0xA7BE68E06D2475BCuL, 0xAFE1A68611CCC795uL)),
                ((-1, -36, 0xA0919859DE85F0BCuL, 0x65A13FAB0FE479C0uL), (+1, -26, 0xC371E31646EDEDF9uL, 0x67DBC45289CF8FF0uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xBC9E2172ABC36676uL, 0x4E4D340B648699CEuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0xE7F8945799CA66ECuL, 0x1A545F44C1A24B00uL), (+1, 0, 0x988B8C709EFEAA63uL, 0x69FA9C3E5E2B01FAuL)),
                ((+1, -5, 0xC89EB0D1EFB2DC22uL, 0xA174E13BF9ACA5E2uL), (+1, -1, 0xC1CAC3AF54595DEAuL, 0xCA530BE01342EF02uL)),
                ((+1, -7, 0xEE6107C36F1CF4A8uL, 0x09ABEE8D6D258EF8uL), (+1, -2, 0xA419DDDDF5F76179uL, 0xE5A8F0FC80A94B6FuL)),
                ((+1, -9, 0xD0C448E99420B9BFuL, 0x72864968D994F944uL), (+1, -4, 0xCAC141F5DBCF28C5uL, 0xE71376BBD5BB2F5FuL)),
                ((+1, -11, 0x8A9A7AEFB2419264uL, 0x178C218101F7F155uL), (+1, -6, 0xBD79E88ABE9FA128uL, 0xBCAF49F8F57B2547uL)),
                ((+1, -14, 0x87DC56360371E8B4uL, 0x6483C42F3D074E29uL), (+1, -8, 0x8765272580C645E7uL, 0x654C89D17353AECCuL)),
                ((+1, -18, 0xB6272974D7C76A31uL, 0xCD261BD74895381EuL), (+1, -11, 0x924E65F1341DCD51uL, 0x2E348074E82A32EAuL)),
                ((+1, -22, 0x9FFB577DD57E1AA6uL, 0x4BB98110662C77C6uL), (+1, -15, 0xE5FB7955D62C4444uL, 0x6A3CEAA3A07A2863uL)),
                ((+1, -28, 0x882F9FD72A92737CuL, 0xFEA63FA0AD071C4EuL), (+1, -19, 0xF027FC7CC904C269uL, 0x23BE86A6DDADAB16uL)),
                ((-1, -38, 0x89949DADB96CEE55uL, 0x8AA95925D746ADD3uL), (+1, -23, 0x815E0F7EEA5B982BuL, 0x3ECE7BE3004353EBuL)),
                ((-1, -41, 0xA107B69F11B61258uL, 0x595F59D39586E49CuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xDBF964A43C64AFE4uL, 0xC81FD7066A2DDAEFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x9C408A5BBBF08819uL, 0xA9C74933AF16C236uL), (+1, 0, 0x9A8A07C397B49AD8uL, 0x341D29423AAF3EAFuL)),
                ((+1, -5, 0x99F1FE13123DC05EuL, 0xBF4145F56FC22C37uL), (+1, -1, 0xC8134AFBA2127F2DuL, 0xAF5C12289FEC6A3EuL)),
                ((+1, -7, 0xC75B0422845C760BuL, 0xA485DD06CA404821uL), (+1, -2, 0xAE551C2110E85155uL, 0x1DD5FE4E90CE85A2uL)),
                ((+1, -9, 0xC8C760C3A8A97B4CuL, 0x567795E2AA13CDCAuL), (+1, -4, 0xE0D0B3078A903470uL, 0x081FD0CEA275B171uL)),
                ((+1, -11, 0x983C15C52CDBE458uL, 0xC118ADE5114ED66CuL), (+1, -6, 0xDF79A385C2719076uL, 0xB22D795E24B7F816uL)),
                ((+1, -14, 0xB591C1610F2C261DuL, 0x8EAA785384E03D85uL), (+1, -8, 0xAE96D21AB38E7926uL, 0x81D9BE333C6CA958uL)),
                ((+1, -17, 0xA2AD2480B5F131E0uL, 0xCE6DD0A85A457FCCuL), (+1, -11, 0xD6C641BCB02F10CFuL, 0x76E79C50B142AB1EuL)),
                ((+1, -21, 0xD78B1CC79A55EB4DuL, 0x6459F4F4807FD63BuL), (+1, -14, 0xCD8C762FB273B516uL, 0xFBA5C31B2FA02684uL)),
                ((+1, -25, 0xBCCF879EE9BD1410uL, 0x1968074374B1D883uL), (+1, -17, 0x9437B256E952512EuL, 0x70BCE67775070A43uL)),
                ((+1, -30, 0xB39402C70FEE2588uL, 0xF056AE45B878E20BuL), (+1, -21, 0x9738C71E2536812CuL, 0x91457B8B5B3093B6uL)),
                ((+1, -37, 0x9BA560D6796E888BuL, 0x899E5571287E904AuL), (+1, -26, 0xBE6D7D7F96F358D1uL, 0xC920422C4BB6BA99uL)),
                ((-1, -46, 0xE0E6B19FFCF92905uL, 0xEF41E9E5288EC418uL), (+1, -32, 0xC8396075A59775C0uL, 0x821F84F521863815uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xC0CE97241CD63229uL, 0xE33FE173B0D3EB30uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xC38B1139179A9DF8uL, 0xA32A7A9124EAAAD2uL), (+1, 0, 0xAD09F896D8D552F0uL, 0x26EE82599809EDB3uL)),
                ((+1, -6, 0xC85712264A32ECC7uL, 0x32717F85A54E359BuL), (+1, -1, 0xE32B958B9D4A9E99uL, 0x7A083D70C80DA245uL)),
                ((+1, -7, 0x841558A4A9C43CB6uL, 0xE1C872A32F080600uL), (+1, -2, 0xBDFE9513F1FFACF2uL, 0x31F9C31C1BB75927uL)),
                ((+1, -10, 0xF5E9B435C0E59625uL, 0x72327DCA88E215D5uL), (+1, -4, 0xE075839442D4402BuL, 0x05DF507EF3CA010EuL)),
                ((+1, -12, 0xA7FB3C0558788E7EuL, 0x2952CCE3148A6138uL), (+1, -6, 0xC4B1A6CC1FA88227uL, 0xD0D3A75AD02AB7D2uL)),
                ((+1, -15, 0xAA904063E2C59B06uL, 0xFAC22603E2073F3CuL), (+1, -8, 0x82D253A9E201FD34uL, 0xCC0A4550DC3226E6uL)),
                ((+1, -19, 0xFF7A80600DACA91DuL, 0xC1A841FB7B4A0EB5uL), (+1, -11, 0x84EA6912F3102FADuL, 0x843DDC7DB186BC01uL)),
                ((+1, -22, 0x8935566AD5E4D2EFuL, 0x1D692581A032F10FuL), (+1, -15, 0xCCB81E112B415134uL, 0xE4D7306757B0C127uL)),
                ((+1, -27, 0xC769B50A49F243EEuL, 0xA53B33B15BCF4E36uL), (+1, -19, 0xE99110DAB3634668uL, 0x9DF3206114C52B15uL)),
                ((+1, -32, 0xAF3610F1F4C19174uL, 0xDA39A40B659CB6EBuL), (+1, -23, 0xBCD58DE9D9F8582BuL, 0x4586077DDE35497AuL)),
                ((+1, -38, 0x92D371B7729D6F5EuL, 0x985BF5DAAE1EBE6CuL), (+1, -28, 0xC79AC646FCA1D85FuL, 0xC1E09D84893F826CuL)),
                ((+1, -47, 0xEA20857D4C0E2991uL, 0xB1053AB7E5FCCA43uL), (+1, -34, 0xEBFB9FB06086FA97uL, 0x4859F87742B3C608uL)),
                ((-1, -56, 0xE30BE839DE97178EuL, 0x174F144E18DD84E8uL), (+1, -41, 0xD8A20934C353BC91uL, 0x3ED67E7608B58C15uL)),
                ((+1, -64, 0x92E5B394AEB0D91EuL, 0x2CCC3911131C0F97uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0x8F3C3123702D81C1uL, 0x0D57EDD710B1169CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -7, 0xDC50243BCDE14397uL, 0x2279CF45E9EB242CuL), (+1, -1, 0xF419B517DB6CBA00uL, 0xED6408FDD5343432uL)),
                ((+1, -8, 0x9A2C249AE8AF8FA1uL, 0x428116079754EA5AuL), (+1, -2, 0xD57DA860A6B3BD3BuL, 0x8EE14232BC053F59uL)),
                ((+1, -10, 0x801D4F00112D9FC9uL, 0x7E8E9238EA43C419uL), (+1, -4, 0xE153DA362753EF13uL, 0x7E21B3351379BB78uL)),
                ((+1, -13, 0x8AD3640EEF690176uL, 0x1B27782E39835801uL), (+1, -6, 0x9EEBD540DED71FDFuL, 0xE8B75CE87110ECB7uL)),
                ((+1, -17, 0xCB23B9C69D43E1E6uL, 0x930288F079C423BFuL), (+1, -9, 0x9CC64C47F62A2749uL, 0xCD14FD3B31364267uL)),
                ((+1, -21, 0xC94CF006833DFC59uL, 0x37B412E37B7D16A2uL), (+1, -13, 0xDC0CF55EA445062CuL, 0x8FC55B3C2991494FuL)),
                ((+1, -25, 0x842D24F15401984FuL, 0x722B7A089071D768uL), (+1, -17, 0xDB405711EFBED588uL, 0xDD7AECD18FE3B3C4uL)),
                ((+1, -31, 0xDB7E13A821B26902uL, 0x04C9D29161F6DE49uL), (+1, -21, 0x981BF35F08ADFA6EuL, 0xD1D46093EB3DE338uL)),
                ((+1, -37, 0xD5201B36D263026AuL, 0x1BC80474F7A18B01uL), (+1, -26, 0x8D71FC784C5074E7uL, 0xFD184FC4C72C7017uL)),
                ((+1, -44, 0xD4C3630B6087CCCDuL, 0xD9C705CA3A1F98EFuL), (+1, -32, 0xA5864AFE23447ACAuL, 0x5764BD4482B84CE0uL)),
                ((+1, -52, 0xAB39918BEE2DA05CuL, 0xAF25EE7E957964C2uL), (+1, -39, 0xDC4AD9D699440CE9uL, 0x303AB1AF64D33C8FuL)),
                ((+1, -63, 0xDA303AD14C0922E5uL, 0xA10146456F242D6CuL), (+1, -46, 0x8BF501D03ADF2561uL, 0x1C42ED691275C127uL)),
                ((-1, -75, 0xF58D1E39164FDE3DuL, 0x51F92F745E55FF32uL), (+1, -56, 0xE7D290B964CEFCDCuL, 0x2D2DD865BCF83E2DuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -8, 0xCC0A874B56266400uL, 0x9BCBEF9585E8B1A4uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -9, 0x99C5A2CD979B3CA9uL, 0x1095E4AF94BDC095uL), (+1, -2, 0xF0D3FB9E87FA1CDBuL, 0xA1EFA24C7439C4D6uL)),
                ((+1, -12, 0xC990E0EB470C9428uL, 0xC5F556670AB218DEuL), (+1, -4, 0xC9A7ABDA1C48F375uL, 0x5202D6288F31B4E4uL)),
                ((+1, -15, 0x96367541B9593443uL, 0x028B983212EE26F1uL), (+1, -7, 0xC5B05BD5CE566E03uL, 0x7FD3A416E272E5B4uL)),
                ((+1, -19, 0x8BCBE981193A2717uL, 0x2E853E89D7BB32AEuL), (+1, -11, 0xFB18ECCABC890BDAuL, 0x4308CDA6BA309511uL)),
                ((+1, -24, 0xA829E238DC0DFEE4uL, 0x72EF313CFF00B112uL), (+1, -15, 0xD7F9CE7C50CA0900uL, 0x8599FBD2BDEFCED3uL)),
                ((+1, -29, 0x82FE0A52B00D9111uL, 0x1C692C09DD911AC2uL), (+1, -20, 0xFF819038FF88157FuL, 0x7670149238AAF766uL)),
                ((+1, -35, 0x810741E976D360CFuL, 0xA5114ABBF27FE6ACuL), (+1, -25, 0xCEEEF43CB46D7776uL, 0x75AB3F72B724E7E0uL)),
                ((+1, -42, 0x9855E28A20102B8CuL, 0x18844EEF873CB6A6uL), (+1, -31, 0xE04466D6785233CEuL, 0x9B5D62E336B72027uL)),
                ((+1, -50, 0xC2F1B396A78F9469uL, 0x07FF40158AD5F531uL), (+1, -37, 0x9B76D1E1E0123F3DuL, 0xBA0E68236DB86975uL)),
                ((+1, -59, 0xDC18DA743658C304uL, 0xE8732A0823AD43ADuL), (+1, -45, 0xFEDF5352351FA8E6uL, 0xABFCD2D7A6F59C7CuL)),
                ((+1, -70, 0xE5F2065F86AC9127uL, 0x1519242F22EF596FuL), (+1, -53, 0xD5CD51C64B3E47D6uL, 0x476F3E9F02E4ED40uL)),
                ((-1, -81, 0x94F21DA5A6C0ABDCuL, 0xA2FDB91E74558E7EuL), (+1, -62, 0x8373C4D69CBB9939uL, 0x274DA795332398B9uL)),
                ((+1, -92, 0x8711D98B60C9FBA4uL, 0xFE3C4F9B0001C18CuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -9, 0x9069CB57E770F408uL, 0x916308F12F0D5C7FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -12, 0xBCE3836CB6FBC58AuL, 0x39A6167DE84D825DuL), (+1, -3, 0xD7686628964B8083uL, 0xCE25B7174B3A42C1uL)),
                ((+1, -16, 0xD39E2175DC5D4E63uL, 0x88D4AFB890D001C9uL), (+1, -6, 0x9F8CFD8B3D08A1FFuL, 0x58AC6AE741E76D35uL)),
                ((+1, -20, 0x84BA7FA739D08F8FuL, 0x3B8309E9C7106A61uL), (+1, -10, 0x88BDBB3E96AB5DE7uL, 0xDFDF92B1541EF01AuL)),
                ((+1, -26, 0xCCABDA87BE52A86FuL, 0x3DCF52EC18CDEECBuL), (+1, -15, 0x95EFA6ECEC263556uL, 0x1D96A91A42C6CBC6uL)),
                ((+1, -32, 0xC8A03B7D406D72AAuL, 0xD510F3D12BF1741EuL), (+1, -21, 0xDBA05CB70A7DC621uL, 0xE02EBA278026F365uL)),
                ((+1, -39, 0xFA2399225517511DuL, 0x6445D1762E068B6FuL), (+1, -27, 0xD9E67A598732F8BAuL, 0xB15B5964D487025EuL)),
                ((+1, -46, 0xC11237C3FE55C351uL, 0x72113A4E12355A0AuL), (+1, -33, 0x9176362217688BADuL, 0x60B950BBF846B67CuL)),
                ((+1, -54, 0xADE3633DB72D8747uL, 0x4757EE98CF8AECBFuL), (+1, -41, 0xFEA13B05D0A752A7uL, 0x437F9340F1FA9366uL)),
                ((+1, -63, 0xA3702C8E5340C071uL, 0xB756023C60C9FE03uL), (+1, -48, 0x8AFAEFE9A127275CuL, 0x48FE418971B37BC2uL)),
                ((+1, -74, 0xFF31F6620B664BE4uL, 0x6D3507CCC5C06B3BuL), (+1, -57, 0xAD6F18EEBFB10C46uL, 0x41BBB62D28E6BB61uL)),
                ((+1, -86, 0xA4BC75EBA927C2F0uL, 0xA6BD9692AF77545AuL), (+1, -67, 0xD2C4BBD6D1B20F6CuL, 0xDECE135CFB1A254BuL)),
                ((-1, -100, 0xC6F7BB39E37256A4uL, 0x30F8BA11C30D09ECuL), (+1, -78, 0xACCDB9A8B0022F15uL, 0xF3380209CB27F2C1uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xCC42299EA1B28468uL, 0x7E59E2805D5C717FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -12, 0xE98327021A46F7B4uL, 0x631C0FFF0AF8AC9BuL), (+1, -10, 0x92550744DFC376E2uL, 0x1FEFAF793553C675uL)),
                ((+1, 5, 0xC42B130CF0BF40F4uL, 0x4A4F450A05229B32uL), (+1, 6, 0xFE9C46A343E32BEBuL, 0xDDD9B88646BF6961uL)),
                ((+1, -6, 0xC547D7AE6B46DF82uL, 0xAE85FFA033D4E352uL), (+1, -4, 0x85A16702E91EC6F0uL, 0x36136F4747F79C1FuL)),
                ((+1, 10, 0xA1F1928511BDD0D0uL, 0xBFF79633E5DFA499uL), (+1, 11, 0xE9609A34A6335532uL, 0x6C5F4117D7C18056uL)),
                ((+1, -3, 0x84B3FF0383B4E251uL, 0x09A273B4235FBA40uL), (+1, -1, 0x881F6CB1FB14BE9AuL, 0xEA8E4099A7E40A41uL)),
                ((+1, 12, 0x89CED50B0594932EuL, 0x3FBF30132191CAD7uL), (+1, 14, 0x9FD7F07652E23C44uL, 0xBF80B20BEC2A1D46uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_0 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD8B2C831775D140DuL, 0x968875EC7736ECC7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x9EAF794633DA12B4uL, 0x217BFB040583C994uL), (+1, -1, 0xCF7845F7CA951E07uL, 0x1C36C6740CC09609uL)),
                ((+1, -2, 0xDA775D29E00384A7uL, 0xF5ED6544847D7084uL), (+1, -1, 0x82B643002E6AE679uL, 0xD77CA1A091168B94uL)),
                ((+1, -3, 0xC42B4A3B4BAB1A64uL, 0x32D0E486F53B9AE1uL), (+1, -3, 0xD50083A84AB517CAuL, 0xF7AEFE0A74DB7BBBuL)),
                ((+1, -4, 0x85769619C7E41504uL, 0xD1BA8E9EDE116897uL), (+1, -4, 0x91DC847B66B8B727uL, 0x764EA67BE1490244uL)),
                ((+1, -6, 0x91DC35C06994AD45uL, 0x4AE4D52D45D9EAAEuL), (+1, -6, 0x9705A86B6F188174uL, 0xA3E15A27D60DF17AuL)),
                ((+1, -8, 0x80970030448F79F4uL, 0xDEDCC075402E90EBuL), (+1, -8, 0x8661950F02992B9DuL, 0x2E7BFCAA2BFDA380uL)),
                ((+1, -11, 0xB6309B6606B0A72FuL, 0xED9F9560B36CED28uL), (+1, -11, 0xB7B7CCD7BC2D1504uL, 0x325A8324CFD54EC4uL)),
                ((+1, -14, 0xCEBC29A3571C20F9uL, 0x28033F6DFABCB7DDuL), (+1, -14, 0xD6382993C2ECCB40uL, 0xC7E8E420C3591C5EuL)),
                ((+1, -17, 0xB56E8C489262E911uL, 0x16E24BCE87444711uL), (+1, -17, 0xB09A8D0E9E945079uL, 0xC2DC213AF2346DB9uL)),
                ((+1, -21, 0xDDB13E9B0E231D9FuL, 0x268149E05B38499BuL), (+1, -21, 0xE8CFDF4102B9A3A5uL, 0x49935EE88F97092CuL)),
                ((+1, -25, 0x8DF96C00F5DFEFD6uL, 0xE4861CF62DDDC83DuL), (+1, -25, 0x829519B04167E730uL, 0x8D3590052E38E727uL)),
                ((-1, -34, 0xAE1F09F850943BB0uL, 0x55695457BA4CBFF6uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA680572CAF1E9B8DuL, 0x6CB90A64B6886A66uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x94E3F1D6AD416BDCuL, 0x7A42E69C824EA721uL), (+1, -2, 0xF0503268DBE9A918uL, 0x03F7359C844D4D61uL)),
                ((+1, -3, 0xE9D9D6B53A893C32uL, 0xD2DC2C61F4195DCBuL), (+1, -2, 0xAEBE5C5CB555BF77uL, 0x6434C5349983D463uL)),
                ((+1, -4, 0xD8D008CB557B74F0uL, 0x75DF45CEF29B6A67uL), (+1, -4, 0xD75AA7A8A9158285uL, 0x25ACBC14506AE661uL)),
                ((+1, -5, 0x8C00ECC7A257C5C1uL, 0xC15A196979D780ACuL), (+1, -5, 0xA97021C179ABFB34uL, 0x2C7D4FB0D7EC3B23uL)),
                ((+1, -7, 0x93047C565F0BB0E7uL, 0xCE6549B720CC8D60uL), (+1, -7, 0x92785796C8A1BB2BuL, 0x10C64E22B0F91A47uL)),
                ((+1, -9, 0x879A5011933C0AF1uL, 0xD2E5D06C023EBBEFuL), (+1, -9, 0x99FC139B4A1988BDuL, 0xB7C3CE3C1E50D15BuL)),
                ((+1, -12, 0xCC65DB4F556D4934uL, 0xF0A83D717D5A8EA9uL), (+1, -12, 0xB9229828C7836002uL, 0x0557DB85C06896BAuL)),
                ((+1, -15, 0xE7DD00752F669639uL, 0xA22B46467E1B2FE8uL), (+1, -14, 0x846033D965F015ABuL, 0x3A123ED7FCB8E9E4uL)),
                ((+1, -18, 0xD45D14E042A42EEAuL, 0x723A6E892B4F3CA6uL), (+1, -18, 0xC9600CF2D3537C97uL, 0x9FE4EAA4250F2E38uL)),
                ((+1, -21, 0xA1D12E610E871D81uL, 0x6B1299910D753008uL), (+1, -21, 0x97C931BE345D4AC9uL, 0xA852C9CDCD7B7954uL)),
                ((+1, -26, 0x9464C26764F995E4uL, 0x23C5F100D5D002B4uL), (+1, -26, 0xCAB2F5B34551ADBEuL, 0xBC5F77727B26F262uL)),
                ((-1, -29, 0xDEC9C24E8EDDA546uL, 0x01A4792F9FAB5700uL), (-1, -29, 0xF38DF28DFDE31354uL, 0x7E35240231143250uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x9692B6FA15E76510uL, 0x6CAA649ED6FDC58CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xD59D79581ED84AEEuL, 0xBC8E4693BF260650uL), (+1, 0, 0xACA8641B7CC1CA29uL, 0xE5C30BF06B2CF208uL)),
                ((+1, -1, 0x8CAA9624553A111FuL, 0x0F0044E12DE8951BuL), (+1, -1, 0xF1B230610C2D663BuL, 0x0055AA7696D30A17uL)),
                ((+1, -3, 0xF4AE1BFA3F9E1BB7uL, 0xEA78C7E74D0A9D71uL), (+1, -2, 0xE0C093B7FFB72687uL, 0xFB502F6989F324CBuL)),
                ((+1, -4, 0x9B13B84C7E74CEA7uL, 0xAF0D68FEFFCF368BuL), (+1, -3, 0x9915C3ECE016911BuL, 0x2ED3331A8965A936uL)),
                ((+1, -6, 0x95E86AA26A59475EuL, 0x6191AF1431E3D53AuL), (+1, -5, 0x9FDD3D1AE4514B4BuL, 0xA6CAD792225C26B1uL)),
                ((+1, -9, 0xE24CBF1B6FC984FFuL, 0x50B03225275CD398uL), (+1, -7, 0x82E4AF4AE1BFC046uL, 0x53CA395355D0A5ADuL)),
                ((+1, -11, 0x85AE3C84913285DAuL, 0x81FF87D41C500CBFuL), (+1, -10, 0xA9324D2B48B0B38DuL, 0x15EF5952BEA5CC41uL)),
                ((+1, -15, 0xF4F43DE86CC08721uL, 0x42A71842D50B3CBCuL), (+1, -13, 0xABB16AC29383EFDBuL, 0x3657007B67BEE100uL)),
                ((+1, -18, 0xA9648ED8F606B98AuL, 0x0BA5B9D2725F2373uL), (+1, -16, 0x8606425BBD42E491uL, 0x59C086585A1FFEF1uL)),
                ((+1, -22, 0xA6B15EA3AFDCA9BDuL, 0x449BF64F7A94B5A0uL), (+1, -20, 0x9A97B63B9533780AuL, 0x91A023523A652022uL)),
                ((+1, -27, 0xD2740857376F5459uL, 0x0C107F330EFDBB10uL), (+1, -25, 0xF3281B1E5261CF64uL, 0x88EE67F71F7B71C0uL)),
                ((+1, -33, 0xF56730B9C0C70E49uL, 0xB77B7C017B2D310DuL), (+1, -30, 0xDBC431AA9DCED510uL, 0xA26293151C411E7EuL)),
                ((+1, -41, 0xB46328B08678ED14uL, 0x69AACEA42B601314uL), (+1, -37, 0xFD5A29B934D9E18AuL, 0x1126DCFADEB8FCC5uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8D2D5192B6EA3F54uL, 0xD42C605DACF4FB7FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xB55C33B6B8E35B92uL, 0x53AC6BE6278084ACuL), (+1, 0, 0xADD38E5BA18EA3FEuL, 0x5E1D3B124DA79276uL)),
                ((+1, -2, 0xE0D2A623348E7D3BuL, 0x76D1B95D28F3F2EFuL), (+1, -1, 0xE5E3F0CDA19A939DuL, 0x2E4C3FE9A0288A69uL)),
                ((+1, -3, 0xB13B9B3E6AAABE2BuL, 0x8D050533BAA40A70uL), (+1, -2, 0xC1D9B74565A7F2B7uL, 0xDFC0FF1A3D4DAE24uL)),
                ((+1, -5, 0xC49315F877F92C5FuL, 0xD5715CAD5B27E629uL), (+1, -4, 0xE6D29AE2BAA1AC4DuL, 0xC4B9BD3A094DBFB8uL)),
                ((+1, -7, 0xA0A805F53E80E267uL, 0xD53FD5519C02C830uL), (+1, -6, 0xCB8E673ED9A1C035uL, 0x996808CDDF57E93EuL)),
                ((+1, -10, 0xC532ED119CB9F4E3uL, 0x57ABB703B38E910BuL), (+1, -8, 0x87D6E673716887F2uL, 0xFA01E64ADAA4EC1BuL)),
                ((+1, -13, 0xB5DAF5ACB8A0BE0EuL, 0x6ACB9CAD0BCF7C44uL), (+1, -11, 0x89C61548FC5A5A91uL, 0x6B69B7F5CB1A4BC9uL)),
                ((+1, -17, 0xF7CC238D377D6AFDuL, 0xE34B85210F5832B6uL), (+1, -15, 0xD214A2FEADCE2677uL, 0xD2545D97710DD995uL)),
                ((+1, -21, 0xF00C8AB61A47AE40uL, 0xF7B128C018EB9AEAuL), (+1, -19, 0xEA10DD7EEF9B0C65uL, 0x718D8E6AC7FC7A25uL)),
                ((+1, -25, 0x99C5C0671182B952uL, 0x2AD15AE4C83FBE0FuL), (+1, -23, 0xB49339BD5182B6C5uL, 0xD9E6F85B2125A402uL)),
                ((+1, -31, 0xE3088CF464A1D23AuL, 0xF64BBF0C6E879FD0uL), (+1, -28, 0xAEC732183CB2E89DuL, 0x1DDCBD38E0F454B5uL)),
                ((+1, -37, 0x92D06FEF2A23D94DuL, 0xDCF1215E4F99EED7uL), (+1, -34, 0xAF8DDBBFC8BFEB9FuL, 0x48572FDA036A04EFuL)),
                ((+1, -46, 0x9D211F73E840A90DuL, 0x859672B084E318F5uL), (+1, -42, 0xF20D090808B93379uL, 0xD1DFD640943E0F75uL)),
                ((-1, -56, 0x97616B8E6E39F238uL, 0x2C4918EB71F68FBEuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD60CB65E5CE493D2uL, 0x3EA3D9F8346A226BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xC1B543174BE71BBAuL, 0x9173C74782BCA617uL), (+1, -1, 0xF68C300976F2C68EuL, 0x86F18AB094698EC2uL)),
                ((+1, -3, 0xA0F6F6A060ADEC18uL, 0x6B6F34E7626E7841uL), (+1, -2, 0xDAD66DAC90C93960uL, 0xEAA54BBBE0E4D494uL)),
                ((+1, -5, 0xA14F150446010C12uL, 0x4A0B64FDF4D4E79FuL), (+1, -4, 0xEB7A1C422F28871EuL, 0xB3856DBF8111B02EuL)),
                ((+1, -8, 0xD77AEC715D4FC037uL, 0x1459199755D242ECuL), (+1, -6, 0xAA18AB5C6BBE2886uL, 0x3E1DA4BDE51A3B6DuL)),
                ((+1, -11, 0xC848ED4B511DC9D2uL, 0x676C1441AF98028CuL), (+1, -9, 0xACAD1DBD94554E11uL, 0xEA887AE88203EE80uL)),
                ((+1, -14, 0x83566ABCCD7F41F0uL, 0xBC64EB652EFE49AFuL), (+1, -13, 0xFAAF5FB70A5BE416uL, 0x5658C0FAD509A3B1uL)),
                ((+1, -19, 0xF16888DAE1B403B9uL, 0x191BC187DA51A4EDuL), (+1, -16, 0x81E3DAD6433949ACuL, 0x77D3AF24A5DBAEEDuL)),
                ((+1, -23, 0x977DD5C7D2FD1705uL, 0x18D29541EEA919CAuL), (+1, -21, 0xBC9BAD96DC794158uL, 0xB5C9CAD37A5D83A1uL)),
                ((+1, -29, 0xF77B69618F8B30CCuL, 0x95BEA4D2B4B3A656uL), (+1, -26, 0xB8D5AF78238A9587uL, 0x25F47AB9B956F46CuL)),
                ((+1, -35, 0xF3606EEAD7360F07uL, 0xA8B4A633D6B46D4AuL), (+1, -32, 0xE5D5E05153407642uL, 0x2FED75282535AB33uL)),
                ((+1, -42, 0xFD7CBACB0F9E9A7DuL, 0xA1C961E1E4D3CF0BuL), (+1, -38, 0xA429CEAB5A4D747CuL, 0x53C800B7885ACC79uL)),
                ((+1, -50, 0xDB8D1EB2385C82E4uL, 0xFABA26184AE1784DuL), (+1, -46, 0xE2E32BD8A969C90EuL, 0xCDDF75400C597219uL)),
                ((+1, -60, 0x9BD1BA4A8646387CuL, 0x36D12A7572BA7AF9uL), (+1, -55, 0xD02F01825A4720D7uL, 0x11F708810AC92F47uL)),
                ((-1, -72, 0xCBA9D9CE143519D3uL, 0x354A0E77C997C9D8uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x98F4E09C7D1EA78BuL, 0x5087D0DED3A69F1CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0xD1D80865D4EBD4E8uL, 0x03F3CDD2E111363DuL), (+1, -2, 0xBF7512C9420E6071uL, 0xB717A7066632E230uL)),
                ((+1, -7, 0xF50926F3F7268907uL, 0xEFE6122B5C76A7D2uL), (+1, -5, 0xF6A91DE089844A21uL, 0x4F55AB22C554492AuL)),
                ((+1, -10, 0x9E0F8CAA612EEC37uL, 0x7E8C555116244782uL), (+1, -8, 0xB24B8F83F09F2187uL, 0xB66212B499F50C11uL)),
                ((+1, -15, 0xF4D8A684DB9EFD69uL, 0x40CB3779454A0DC5uL), (+1, -12, 0x9DF618EC6C9D745FuL, 0x29FAE2DCF661A15BuL)),
                ((+1, -20, 0xE8038C55BDCC296FuL, 0x7CCCB5E15B1360D3uL), (+1, -17, 0xB0031CBC2416E092uL, 0x941546A3C862FB2BuL)),
                ((+1, -25, 0x8381D7C3BA634387uL, 0xA0EE6A8355E097ACuL), (+1, -23, 0xF3C4EAD57C00C500uL, 0x1FED8ECFB3FD8C3AuL)),
                ((+1, -32, 0xA7289CB5A2D698C2uL, 0xA371994B943FA79CuL), (+1, -29, 0xC82877EDE798C620uL, 0x6D3FB8CFC86DAF7DuL)),
                ((+1, -40, 0xD11AE12CDBD738D1uL, 0xB2F59716A91B5311uL), (+1, -36, 0xB11A20A051255603uL, 0x64699E5E131DEDF8uL)),
                ((+1, -49, 0xC304BE661827ACC6uL, 0x28D286AEFD40101AuL), (+1, -44, 0x8BDC436EB2480562uL, 0x7104599870EC01C9uL)),
                ((+1, -61, 0xE48C9342FC26E4D0uL, 0x885AB201241BFF0EuL), (+1, -55, 0xFB2BE8A035A3AF33uL, 0x8D6F5A86DA26E975uL)),
            }));

            private static ddouble PlusValue(ddouble x) {
                Debug.Assert(x >= 0);

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
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    ddouble u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u;
                }

                return y;
            }

            private static ddouble MinusValue(ddouble x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 2d) {
                    ddouble y;
                    if (x <= 1d) {
                        Debug.WriteLine("pade minimum segment passed");

                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0);
                    }
                    else {
                        y = ApproxUtil.Pade(2d - x, pade_minus_2_1);
                    }

                    return y;
                }
                else if (x <= 32d) {
                    ddouble v;
                    if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else if (x <= 8d) {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }
                    else if (x <= 16d) {
                        v = ApproxUtil.Pade(x - 8d, pade_minus_8_16);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 16d, pade_minus_16_32);
                    }

                    ddouble y = v * Exp(-(2d * x * x * x) / 27d) / x;

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
                ((+1, 0, 0xD9F45C5BF647CCDFuL, 0x867230D3E9C07E73uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 5, 0xDBA44749B3F8BAB2uL, 0x60143AA44F485171uL), (+1, 5, 0xA4E5F58EA73A74F6uL, 0x3C9991377F6A55A1uL)),
                ((+1, 9, 0xA8756B2D0C968E20uL, 0x3654463900FC64B5uL), (+1, 9, 0xB3A061B32C420AC2uL, 0x36BA79D94C2FD7DFuL)),
                ((+1, 11, 0xD55D3BE4D354A1B1uL, 0x7BCAA2518434961BuL), (+1, 12, 0xD560C48D2F34A1E8uL, 0x96A3B7A2258FC2ACuL)),
                ((-1, 9, 0xE4B4EE934E9DFFE5uL, 0x5E556EAA73A223A2uL), (+1, 15, 0x93787FE3C1B956B0uL, 0x787C2AE315000E98uL)),
                ((-1, 16, 0xAD0FA0D4966857B2uL, 0x5602039105A19A26uL), (+1, 16, 0xE49CEBB8B0475231uL, 0x2CCE1D0C20E6485CuL)),
                ((-1, 18, 0xB7674F0E0D3DCAECuL, 0x4F62D291D6648ADEuL), (+1, 17, 0x994594AFF343D7D1uL, 0xC612F93EC09D114EuL)),
                ((-1, 18, 0xF2CF827446CA210FuL, 0x05F2D85EFE85D538uL), (-1, 16, 0xBBEE2B9995C2435AuL, 0xBE64BF014EE825D1uL)),
                ((+1, 18, 0xCA99B548EBAE162EuL, 0xAB2231D663757CAEuL), (-1, 18, 0xEB8FA60684D56B59uL, 0x228860ACA1AAD6B2uL)),
                ((+1, 20, 0x9F024B6432E29C36uL, 0xEA75324D1930238AuL), (-1, 17, 0xC83E604476E80438uL, 0xB1A24378DAC583F5uL)),
                ((+1, 16, 0xA2E8FDF123DD9AA3uL, 0xBDA4E4A3ABE231FFuL), (+1, 18, 0xB75F4E127491CC76uL, 0xF543F3FABA2D0AE8uL)),
                ((-1, 19, 0xC7B964351AAAF7AEuL, 0xD2CF44B2202C61F2uL), (+1, 17, 0xA299CF9BA11F5384uL, 0xD6791C0B0E5992E0uL)),
                ((-1, 15, 0xBE2E1E761B2B4CD9uL, 0x1AC991315C203AB4uL), (-1, 16, 0xAB35378798C43952uL, 0x5D81D7CB8D6ED4E8uL)),
                ((+1, 16, 0x9D447E48A6601E03uL, 0xD38EBC873544A187uL), (-1, 13, 0xAC7C7E09C2857ADBuL, 0x14BEE8AF047CE4E2uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xF688C14B5900BCE5uL, 0x754AD54726E72685uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0x8136F88D5DAD8622uL, 0xAF0C0D312EB7325FuL), (+1, 3, 0xBCD4CFFAFE11D6EEuL, 0x767441B7E7F269CCuL)),
                ((-1, 5, 0xAEA08A71FB684917uL, 0x4ACC315239F027EFuL), (+1, 5, 0xB0FAD6532443997EuL, 0xA146B588C2E084D5uL)),
                ((-1, 7, 0xA51C394B9EAF8E36uL, 0x6E92D7FF27C839FEuL), (+1, 4, 0xCFEAF6E539451258uL, 0x9483B9EFE0CD82A6uL)),
                ((-1, 4, 0xDF2FACA21C525815uL, 0x5119C0E0F964CF9AuL), (-1, 7, 0xAA53101E52575B50uL, 0x78EAD9F2DFEB2D63uL)),
                ((+1, 9, 0x9F6894CFE5660275uL, 0x6A78883B94F0596CuL), (-1, 7, 0xE4C11D64CD8410B2uL, 0x3EBFED7116798FECuL)),
                ((+1, 8, 0xD4A10EF0F9E1DE02uL, 0x84B2912DE36DBCBFuL), (+1, 7, 0xE545F5FB14A54DFAuL, 0x07086E86D9DC6906uL)),
                ((-1, 9, 0xDB4F8D79DBDE0926uL, 0x6665CA2AA077BB13uL), (+1, 8, 0xAC6043903857A3E9uL, 0x00AA630A0226463EuL)),
                ((-1, 8, 0xE8D1436CBC67B791uL, 0x4FBF86AB57EF3D18uL), (-1, 7, 0x968BDE89858DE17DuL, 0xCA8FDD4570781F96uL)),
                ((+1, 8, 0xE43F677A678648E2uL, 0xE92F4FB4916E6902uL), (-1, 7, 0x92301BAF4CAA3FE7uL, 0x9FD924E1FFFB3E18uL)),
                ((+1, 6, 0xABF4F3A19D4D5DBFuL, 0x4A7E45B32B7DB388uL), (+1, 5, 0x94BDCFC09A6D32E8uL, 0xF9A3216E5FCFB9E1uL)),
                ((-1, 5, 0xA05A7C479647EFE6uL, 0x32A24437FE5489BCuL), (+1, 2, 0xF1102878C99A1AD8uL, 0x54252EDC8DEB95E1uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -3, 0xCD49451AD609DA8BuL, 0x7FA41165FF952667uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 2, 0xB3F128581AD38CC5uL, 0x6244088FEA30B376uL), (+1, 2, 0xA5A7258004406BE0uL, 0xBD45485C6E67FD32uL)),
                ((-1, 4, 0x96E4FFB6640CAED9uL, 0x2F02D2986268745FuL), (-1, -5, 0xA7C1721ADE3D76E5uL, 0x85099946DE9E82F2uL)),
                ((+1, 4, 0x9F0BB8A39F4802FAuL, 0xE878D3E77F579D21uL), (-1, 5, 0x8567623E15E3E561uL, 0x9EF51C4766C9240EuL)),
                ((+1, 6, 0xEF19C26BFA4F8966uL, 0xF68EF019D5C7EC27uL), (-1, 4, 0xB678763EF3EE17B7uL, 0x84CB3E31B4B032F6uL)),
                ((-1, 2, 0xA3B33F018D2D1F40uL, 0xE35B5CD7911FA974uL), (+1, 6, 0x9E3C30D85C2DC03EuL, 0xC4344F2386DF6799uL)),
                ((-1, 7, 0xFFC041426F81E6AEuL, 0x04912624773C2FF1uL), (+1, 5, 0xE86205DDF06458F5uL, 0x551D5DF5187E2096uL)),
                ((-1, 4, 0xE1CB6C46644CF770uL, 0x8D109A49CE3B0DEEuL), (-1, 6, 0xA502BE76AE5B199AuL, 0x80ECCEB394E94448uL)),
                ((+1, 7, 0xD11E52F00A47E2EDuL, 0x1B357E64D18B9C74uL), (-1, 5, 0xAF91F0801E57760BuL, 0xA638BE39A4177BAAuL)),
                ((+1, 3, 0xCE24B57238C66228uL, 0x3C2964CBB8C37618uL), (+1, 5, 0x817DD974FED07EE5uL, 0x821E8E6F44ABD311uL)),
                ((-1, 5, 0xBB1B9E643283AA8BuL, 0xF98E6A7C8E7C489AuL), (+1, 2, 0xF0430E0405A6B74FuL, 0xDC464960F57A2A6DuL)),
                (Zero, (-1, 0, 0xFD354718481A34DEuL, 0xE1073E13EF8C4E5EuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD9F45C5BF647CCDFuL, 0x867230D3EB04A27BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x806625577A4E8125uL, 0x953788A0199861A0uL), (+1, -2, 0x8BFDEE9BF9DDE137uL, 0xAC3F458540D27333uL)),
                ((-1, -7, 0xE3236A7953C01320uL, 0x1C66B22C55945F6FuL), (+1, -5, 0xC5E66882630ABC57uL, 0x28C5F47EA1E18E9BuL)),
                ((+1, -5, 0x80E13E5BEB9FA20DuL, 0xBD04309A38C7E943uL), (+1, -5, 0xD1488AB87256CF33uL, 0xEF0AFED490C39A79uL)),
                ((+1, -7, 0x942913B9F489D427uL, 0x8EDF0E61BA382BA5uL), (+1, -7, 0xAC5F0B0C70FDF498uL, 0x4DBDED0A12ED4E52uL)),
                ((-1, -11, 0xD389EBE17C94EA16uL, 0x5CDDB816503F3C3EuL), (+1, -10, 0xD8D03DA68B05D6BDuL, 0x2B7977525B75DE4AuL)),
                ((+1, -11, 0xC48B7E7F147F51CAuL, 0x77A1C269E46688BFuL), (+1, -11, 0xBFFFAE08AB795076uL, 0x1F9BABA196CD3FD2uL)),
                ((+1, -15, 0xC936734A21976E53uL, 0xA0E268292DD4E2D4uL), (+1, -14, 0xE422D99524B34BFAuL, 0x3066F50BD9B3CCAEuL)),
                ((-1, -18, 0xC36196447FF05F4CuL, 0xEB3973483D7BF6FDuL), (+1, -17, 0xE70C02917B94F030uL, 0xB6FABA095890DAB6uL)),
                ((+1, -18, 0xB0107D12B8D4F63BuL, 0x7ACFE55343BF5118uL), (+1, -19, 0xAED3C2ADC193A2C2uL, 0x578B5D834BB45373uL)),
                ((-1, -22, 0xD6A470D9D5AEDEFCuL, 0x53620B012915AF93uL), (+1, -22, 0x8C5021D4012E9D26uL, 0x4DFC0111426E239EuL)),
                ((+1, -25, 0xD806767F0576BC37uL, 0xF74855D4096E22D7uL), (+1, -26, 0x89401370D7640A6EuL, 0x6F71406DBA8BB6E1uL)),
                ((-1, -30, 0xC16C422322F31883uL, 0x8CE17C652F6B9593uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm4_6 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x81E11B300568573FuL, 0x720A0A7C1D82B984uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xE5DED5F61D206EA1uL, 0x436850AA2ED137AEuL), (+1, -2, 0xB65C8E4D384F8FF1uL, 0x59618DA035FDD9C4uL)),
                ((+1, -3, 0xAAE18DC0B75AADEDuL, 0xEC846F58B171C112uL), (+1, -2, 0xB2EB3F5D95708281uL, 0xF06396A0C424177DuL)),
                ((+1, -4, 0xB62883B2ACB27089uL, 0x844205A221030F64uL), (+1, -3, 0x9489B4F197F73C8BuL, 0xCECD30550F83F271uL)),
                ((+1, -6, 0xF41717A6D84ABB2EuL, 0x79596F38F58CC323uL), (+1, -5, 0xFD40BFD0E4CCC15AuL, 0x24DEFD9725AA37CDuL)),
                ((+1, -7, 0xDBAA67FF0EEB1D56uL, 0x26D465D43A822E4DuL), (+1, -6, 0xBC22E173EEAEBAE7uL, 0xC79B8B21A570041DuL)),
                ((+1, -9, 0xE76741CE079F1A7CuL, 0x40F6125881FA4EF9uL), (+1, -8, 0xDF9700646101F33EuL, 0x9D869FA31B4DA4CAuL)),
                ((+1, -10, 0x8A2851670C25A100uL, 0x1955E0C7EC16E1E0uL), (+1, -10, 0xFCD74C99C829F35BuL, 0x4800A10C81C6A21BuL)),
                ((+1, -12, 0x830BCFF2D46A3EBDuL, 0x63E6154EBBAE1A17uL), (+1, -12, 0xEC06CA276DC63BBBuL, 0x76326C99C409A997uL)),
                ((+1, -15, 0xC7840F362058ACFCuL, 0xCA634C3C9D000558uL), (+1, -14, 0xC21090628B37DC32uL, 0x74E15EC4BABE418FuL)),
                ((+1, -17, 0x9E492781ABADFEF4uL, 0x8907E100D6785382uL), (+1, -16, 0x88482EBD195018C3uL, 0x1CB7EE3E906E7C66uL)),
                ((+1, -20, 0x9A1131381BDD39A0uL, 0x4A0F5F16D0799B31uL), (+1, -19, 0x9BE2C8CE99236F3DuL, 0xB700E7DAF75E5CF7uL)),
                ((+1, -23, 0xAB2899F96D64E12CuL, 0x595EADA0C4E966D7uL), (+1, -22, 0x8F72BDBE9E13EFCDuL, 0x493EA8723AAFFD16uL)),
                ((+1, -27, 0xB2F1E4174234173BuL, 0x3621072CAD45C237uL), (+1, -26, 0xBB7C129E352AA2C4uL, 0x7B4289155916DB48uL)),
                ((+1, -31, 0xB873D78A5E27FE36uL, 0x26E191F27F3BACE4uL), (+1, -30, 0x9214B15FFEB8214CuL, 0xAE791F11825F3AF7uL)),
                ((-1, -38, 0x8C601203C622A803uL, 0xFDF296DB208328CBuL), Zero),
                ((+1, -44, 0xB1E2F2B9AAFED9F8uL, 0xA3B473BE5D150C8BuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8A1E5174B0993855uL, 0xCCAF23E9CC6EF5A4uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xD7B12A5E2B576F0AuL, 0x82FB8E4E92B389A9uL), (+1, 0, 0xC71BC73589A35729uL, 0x4A9D28A3611C9916uL)),
                ((+1, -1, 0xB92243B81019B407uL, 0xC16BB03D70BC523EuL), (+1, 0, 0xAAE2F4A42CD21F54uL, 0x8F1B0CF021BD218AuL)),
                ((+1, -2, 0xD5B5EE92E61A7519uL, 0x951D559D52480A9BuL), (+1, -1, 0xC525B3374CC2136BuL, 0xC918713B1A276AA7uL)),
                ((+1, -3, 0xB58C0B882501C2CEuL, 0x3CDBAFA67CB6C6E8uL), (+1, -2, 0xA7829C400871AA18uL, 0xB87A7FFF9A588B28uL)),
                ((+1, -5, 0xEC430B4F5BE42362uL, 0x072B1580909FC955uL), (+1, -4, 0xD9F609D0F586EDF6uL, 0xE58E32120D593C80uL)),
                ((+1, -7, 0xF08E387F278F91AEuL, 0x89D810E140C06339uL), (+1, -6, 0xDDF8002BFF03B5E7uL, 0xB3AC523EF14BB71BuL)),
                ((+1, -9, 0xC19AC50641C06CBDuL, 0xB489A48029984B9EuL), (+1, -8, 0xB29C7E793F2BE6A5uL, 0x06D964C7DDAA0B19uL)),
                ((+1, -12, 0xF6B2F4CF8DEBD6D8uL, 0x6C5365972B06747FuL), (+1, -11, 0xE3A2DD8EC6D0CE7CuL, 0x907AFCBDA10074EBuL)),
                ((+1, -15, 0xF7C5241EEDAB1601uL, 0x339D7853DE8F4B77uL), (+1, -14, 0xE48FA86306F4F818uL, 0x4B3663BC234EAFD9uL)),
                ((+1, -18, 0xC0DDE60D2B0F1C0FuL, 0x388FCAD1131202E5uL), (+1, -17, 0xB1FF39CC52284A33uL, 0x685652A809617E80uL)),
                ((+1, -22, 0xE14D897A3888B091uL, 0xD146499E7E572D68uL), (+1, -21, 0xCFC98BF64218D07CuL, 0x3C1DA74EBBDD4CE7uL)),
                ((+1, -26, 0xB648A0ABCDD367A7uL, 0x3FA4D1702735704FuL), (+1, -25, 0xA84536B467767FCFuL, 0xB484321C86E9E324uL)),
                ((+1, -31, 0xA40072EAC4FBFFDFuL, 0x71A023F834B804D2uL), (+1, -30, 0x9734654C6F3AE6C5uL, 0x24952DB454933701uL)),
                ((-1, -47, 0x8697A9D9DC64213FuL, 0x59734C7DCA0CF9B2uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm8_12 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8AB1BDECE0E35810uL, 0x1E332837BCB2CF08uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xD9E87428BB0BA2F2uL, 0x34B0FE44F28C77F6uL), (+1, 0, 0xC90E64080CE11A69uL, 0x63B9AE82CE152738uL)),
                ((+1, -1, 0x9D7085F6D1F5E3B4uL, 0xD5E28AC97A3EC419uL), (+1, 0, 0x91419D1C416ED799uL, 0xFDA43FC488F601B8uL)),
                ((+1, -2, 0x8C455BFE1D5BD9A3uL, 0x415818B3383FB87CuL), (+1, -1, 0x816B318287167CCCuL, 0x54DE4F3A69E8335FuL)),
                ((+1, -4, 0xAD925FA2C1161C39uL, 0x8BD34A22DAEF63A9uL), (+1, -3, 0xA024E3ED93C5DD06uL, 0x123B5BE0AD8938FAuL)),
                ((+1, -6, 0x9F36A7267A70DDF6uL, 0x11D5F864E0D5F747uL), (+1, -5, 0x92E55BE8B83FFD75uL, 0xFD433767D325CAB4uL)),
                ((+1, -9, 0xE17B4696C8494D8BuL, 0x012F0E1E62838999uL), (+1, -8, 0xD0094F915D4CF914uL, 0x3D11DE9DCB1B5119uL)),
                ((+1, -12, 0xFCFE090833E655EFuL, 0x6BA95F8BCDB2A0E9uL), (+1, -11, 0xE96BBD49070C7A00uL, 0xAD78456CDD7926A6uL)),
                ((+1, -15, 0xE469916DB376911AuL, 0x63DD18825C8D868DuL), (+1, -14, 0xD2BD812E2F8B6DAAuL, 0x57A4C7B9273D4D1BuL)),
                ((+1, -18, 0xA711A720DA324D52uL, 0xCD64DB8351CFFAF3uL), (+1, -17, 0x9A24E8B9849C932DuL, 0x339C1939E787D138uL)),
                ((+1, -22, 0xC58A1707102D5817uL, 0xF5DA1873DAE04EA4uL), (+1, -21, 0xB641AA8B4A8C08E8uL, 0xE5C52539CBD6807AuL)),
                ((+1, -26, 0xBA18101D79E39917uL, 0x429D8CDB5222B571uL), (+1, -25, 0xABB26859E6A32200uL, 0x4944E61B8B92FE82uL)),
                ((+1, -30, 0x871B9EC099881908uL, 0x8F22923DFF1BAFF1uL), (+1, -30, 0xF94F85B18CB30A5CuL, 0x63E55A5F5D341C2EuL)),
                ((+1, -35, 0x8BD52A41C91CB180uL, 0x60922361B4C8CAAFuL), (+1, -34, 0x8103B47565DC1EB8uL, 0xE42CA87BDED422F1uL)),
                ((+1, -41, 0xAA9E572BBFDCDD7FuL, 0x10B7B5E32CD2F3EEuL), (+1, -40, 0x9D6B2EB1392E790AuL, 0x29C1FE55022D8A84uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm12_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8ABBA0CDB01EF751uL, 0xA8175C89F86A7D76uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xC5409D5634507A97uL, 0x99FB8947458F9555uL), (+1, -1, 0xB5FDE63FF67F4757uL, 0x489C8E7FDD102BA0uL)),
                ((+1, -3, 0x8634EB7C25178469uL, 0x9D94DBA0424C794EuL), (+1, -3, 0xF7A5DC5ED72D7E8BuL, 0x9A49FBBC41B916AEuL)),
                ((+1, -6, 0xE77440624536FA7DuL, 0x84FD0C118FEC79B3uL), (+1, -5, 0xD58C27AF1FD36FC0uL, 0xCACF9F3D30361BBDuL)),
                ((+1, -8, 0x8D0D3D820651A412uL, 0x9FF64B13DCB9145FuL), (+1, -7, 0x8223AE2B289B8FEEuL, 0xB35F92694B78EE26uL)),
                ((+1, -11, 0x8020B17EEE106F5FuL, 0x382579823184A264uL), (+1, -11, 0xEC6E17738775036DuL, 0x1F7F6042E3BFE624uL)),
                ((+1, -15, 0xB1F09571A6A562FAuL, 0x627E05297A615F74uL), (+1, -14, 0xA42C715F88732D1FuL, 0xDD1E87A4EEC86DBAuL)),
                ((+1, -19, 0xBE1602E77832E4AEuL, 0x664968B03A8DF613uL), (+1, -18, 0xAF613F56545531E6uL, 0xC68C7C02B0FA1901uL)),
                ((+1, -23, 0x9A466EFC867B195BuL, 0x444373647E38CC96uL), (+1, -22, 0x8E56F91CC0B6B51DuL, 0xFF92F1721D92CD98uL)),
                ((+1, -28, 0xB70C8B6E1632A8B5uL, 0x8136E5AB53D18A7FuL), (+1, -27, 0xA8E327D70DD3D37AuL, 0x8860ABCE3E92C972uL)),
                ((+1, -33, 0x90E8DE0F47842D3CuL, 0x45909487748DA964uL), (+1, -32, 0x85B2E466D4C57D89uL, 0x96C372A335B80297uL)),
                ((+1, -40, 0xED3A4FF6852C1069uL, 0x4F4BF16090DF2732uL), (+1, -39, 0xDADFDF83DFCE0A9CuL, 0x2D8693FFCCAB8331uL)),
                ((-1, -67, 0xBC7DD0FFEF7FAD0BuL, 0xA1CA216B2F736AD4uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8ABBAAB22AAF1EA9uL, 0x330E62572555E612uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xDF6F1F53D55EF91DuL, 0x55EAAA3A9A99BC43uL), (+1, -1, 0xCE25EA9684F97483uL, 0xDFB1CDC1CC9898B9uL)),
                ((+1, -3, 0xAFE88213CD3D6348uL, 0x116D2A1E5088BE8EuL), (+1, -2, 0xA24C9251DCCA6814uL, 0xD1949E256ACEBA60uL)),
                ((+1, -5, 0xB4212E816639C2E3uL, 0x06B6D3C426350EACuL), (+1, -4, 0xA631A3838D56C533uL, 0xAFADF13FE3E5BA9EuL)),
                ((+1, -7, 0x869817CF9BDED9BAuL, 0xC32F2132A78210CDuL), (+1, -7, 0xF85CCCC8B26AE0EFuL, 0x01AD469E1F7C7996uL)),
                ((+1, -10, 0x9C1135A5D2D59DBBuL, 0x28652C032C3C7365uL), (+1, -9, 0x8FFE3B81C50AD459uL, 0xA766BEFCDFB94A10uL)),
                ((+1, -13, 0x91B794AC3C7678CAuL, 0x0537C146B2D48E52uL), (+1, -12, 0x8671978E4B92F939uL, 0x0495CE32FC3C87DFuL)),
                ((+1, -17, 0xE055B91FEA65D018uL, 0xB08234127294B5ABuL), (+1, -16, 0xCEFAAD4728744BF0uL, 0x0453F944C70A5E2AuL)),
                ((+1, -20, 0x9081EEEBEFA0B9A1uL, 0xA8E4D5DFF1E489A0uL), (+1, -19, 0x8553E6802E444612uL, 0xFDC903917F9A4075uL)),
                ((+1, -24, 0x9D200B57870B871CuL, 0xB25432FF6B090B8DuL), (+1, -23, 0x90F81D35DCB293E4uL, 0xA2DE524A4755B4CCuL)),
                ((+1, -28, 0x909DDC5204F4C5F7uL, 0x5958B9FDAD21A2C1uL), (+1, -27, 0x856DAAC92D4DE4FCuL, 0x879F13354261223AuL)),
                ((+1, -33, 0xE0B8CC04A56E1037uL, 0x81DD7A816C33FC51uL), (+1, -32, 0xCF5615F9EA0307C0uL, 0x2DBACA82718718E5uL)),
                ((+1, -37, 0x91F3417B196F537EuL, 0xD88A3A8E5BFE4148uL), (+1, -36, 0x86A8A678A01BE11AuL, 0xB5B5F04978AA6584uL)),
                ((+1, -42, 0x9B48A3D01FF633FFuL, 0x4D9E9FA0AEE63617uL), (+1, -41, 0x8F452E0B86B5CCC1uL, 0x6285D084B772EA95uL)),
                ((+1, -47, 0x8222D4F6A4700BDEuL, 0xA029FB12639BE227uL), (+1, -47, 0xF022DD62F122E20AuL, 0xA389F56CE8C26686uL)),
                ((+1, -53, 0x9E49658BFE6B227DuL, 0xCB4B5B7D600F76CAuL), (+1, -52, 0x920A763A37976096uL, 0x73A8A5923D85AE19uL)),
                ((+1, -60, 0xE2259796CCFCB3DFuL, 0x84B777074C16AA3EuL), (+1, -59, 0xD0A6A89FFDEBE518uL, 0x0D824863726ADD76uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8ABBAABC1919B41BuL, 0x619EDF592209438DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xFA64AFA9505CC112uL, 0x4A2720BF50FED2E0uL), (+1, -1, 0xE7058A96422C0DBCuL, 0x9F939595AE90DBFCuL)),
                ((+1, -3, 0xDE55A467A36C9A62uL, 0x9BA8CB2102DC1AE0uL), (+1, -2, 0xCD22368A451BD7F2uL, 0x8401258C901AFE30uL)),
                ((+1, -4, 0x814D1663E886178DuL, 0x2CC0865CC786174BuL), (+1, -4, 0xEE9872D63875A4BBuL, 0x2E5EA6D3A2F5A8EAuL)),
                ((+1, -7, 0xDD3A07CB9CFFCEA6uL, 0x3988EB2032AC45F3uL), (+1, -6, 0xCC1C8AFB199FB644uL, 0x0E3FD03C57EB19B5uL)),
                ((+1, -9, 0x94283CE12022AB10uL, 0xECE9BD6E85C65364uL), (+1, -8, 0x88B1EC29AC08485DuL, 0x9AC5F06D857C8169uL)),
                ((+1, -12, 0xA168CF9121FA8AC0uL, 0xD50FCFBDFC93D6E6uL), (+1, -11, 0x94EC07838EC6AB22uL, 0x3A8192D361B23F32uL)),
                ((+1, -15, 0x92A94CDFE425ABFDuL, 0x3B00FE93DA9CBA90uL), (+1, -14, 0x87509C658FE246E1uL, 0x0384696179B40FDCuL)),
                ((+1, -19, 0xE1F9E8946B22F03BuL, 0x54A273ED83ACA6DAuL), (+1, -18, 0xD07E5AC9D0ED3EF8uL, 0x63D8873A3506710AuL)),
                ((+1, -22, 0x956CACD306C35727uL, 0x3D611A87B388492CuL), (+1, -21, 0x89DD427D9E899534uL, 0xB37D584B59C1F4C3uL)),
                ((+1, -26, 0xA96C2CF0467E8AB8uL, 0x1F6B7A6FB6628D8AuL), (+1, -25, 0x9C50B0BE5220822BuL, 0xFC79893E65738B28uL)),
                ((+1, -30, 0xA9C9B5D05356C668uL, 0x22BF92BDFE5282FAuL), (+1, -29, 0x9CA6FD20099DFB98uL, 0x84E968DEFB9EA744uL)),
                ((+1, -34, 0x86883548B2885416uL, 0xD1094847F749A0E6uL), (+1, -34, 0xF83F7CED3A1ADDD1uL, 0xC717923873C8FFBDuL)),
                ((+1, -39, 0xE6F84AA2D88EABCAuL, 0x8408DEC3CDDD18EBuL), (+1, -38, 0xD519D5EE3A1EE6C8uL, 0x85546C7FD2BF4572uL)),
                ((+1, -44, 0xA28919572C32E9D7uL, 0x297AE87D27434C10uL), (+1, -43, 0x95F6039E672CD26CuL, 0xDAFA823B494A2A55uL)),
                ((+1, -48, 0xB6DD04F1C7653EFCuL, 0x9E0EE355B72AE34FuL), (+1, -47, 0xA8B75563CF986446uL, 0x2C73FCD22EB50EBFuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm64_80 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8ABBAABC1919B425uL, 0x500976719C4496B6uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xB13ECCEB3615BB57uL, 0xD598AB858218F27FuL), (+1, 0, 0xA38861E4DDFC7A11uL, 0x7DDCE0D6C436509CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0x8C027C97F1540D45uL, 0x0EA7EC7DE59F88A6uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 6, 0x8A05C7B34FCFF7B8uL, 0x163497F6CF05A5BFuL), (+1, 5, 0x8806C8D7B0224107uL, 0x0556CEBBC1A22215uL)),
                ((-1, 9, 0xD477B0A412573BF6uL, 0x5FF678617B3269A0uL), (+1, 8, 0xE93B0FB77D7F619CuL, 0x14EFC23A74CD56B4uL)),
                ((-1, 12, 0x9972A09D85821DEEuL, 0xD88D0B37188EC5D2uL), (+1, 11, 0xC9C6132EDF42BC4FuL, 0x290C720FC195E025uL)),
                ((-1, 13, 0xA78B6FD421D98573uL, 0x479C2BA9EF59713CuL), (+1, 13, 0xAC72184659DA5038uL, 0xE1D1DBFE5202A9D8uL)),
                ((+1, 13, 0xF7D3F0CCEF66FB0DuL, 0x0722B8F90F04AEDDuL), (+1, 13, 0xA8CA02E5A101EF9BuL, 0xD2EC4DF6AB157788uL)),
                ((+1, 16, 0xC805C840F1DCD687uL, 0x15C863E07D756F1DuL), (-1, 15, 0x8D854B797A789D27uL, 0xB4B2A5D2A7DD7CD7uL)),
                ((+1, 15, 0xC7502046588CB057uL, 0x514E999AFFEA7F36uL), (-1, 16, 0xA6DCE4CE9BF1671DuL, 0x151EF0B1650656C4uL)),
                ((-1, 18, 0x8487C1AD1C4F8B52uL, 0x4D9A1456D4A32D78uL), (+1, 14, 0xCEA7D92D81E940CDuL, 0xB4B87293BF510B4BuL)),
                ((-1, 17, 0xB2E48CDF3EE9452FuL, 0x70F55FDAB07BDFD7uL), (+1, 17, 0x8E6FAF90E0D09CE1uL, 0x60271E82DDA701C0uL)),
                ((+1, 18, 0x94887A881D8C302BuL, 0x1DF6C2E7DE6E2D45uL), (-1, 13, 0x929BAF108B1216C3uL, 0x5F4A053F0A1A1721uL)),
                ((+1, 16, 0xC719667C02DA6E9EuL, 0x746DF7675A0C37B8uL), (-1, 16, 0x9AD44F09D342280CuL, 0x1C53FB4429CA55EAuL)),
                ((-1, 16, 0xCC253046010C0017uL, 0x685EC3A281258326uL), (+1, 13, 0x93059F1B637E4507uL, 0x46D017E3311C0D22uL)),
                (Zero, (+1, 12, 0xA48B5AA43A8F8229uL, 0x6B4B7012B0CBF647uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 0, 0xD0FFFF3D2094B0F2uL, 0x173B24E2B93482DFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 3, 0xF60C4E35563E41BEuL, 0x9BDCDEDC0F1E3EEBuL), (+1, 3, 0xBCA5585658F2B986uL, 0xB00BF1A530F8B0CAuL)),
                ((-1, 4, 0xEF11D60623859551uL, 0x37730F5FBEA56B03uL), (+1, 5, 0xB2BE1C5AA9EE7AB0uL, 0xD77AD373567FA49AuL)),
                ((+1, 6, 0xD2CF33BBB0AED1CCuL, 0x1ACFDD121DF5F55BuL), (+1, 4, 0xE851079BBAEDB85CuL, 0x88F7B920FEDD318EuL)),
                ((+1, 8, 0xB95FBC64DE4FF508uL, 0xD3194384729969E1uL), (-1, 7, 0xA9F142785BE4C5C0uL, 0xAD4C0F5B4D626228uL)),
                ((-1, 6, 0x91B2A9FA88F57973uL, 0x87DFC72412C3E530uL), (-1, 7, 0xFB9A32EA01DC8165uL, 0xA4E84E36852FA865uL)),
                ((-1, 9, 0xF6793960CBC1B864uL, 0x19721A5BEDE828B1uL), (+1, 7, 0xD540C48B7FD1DC06uL, 0x0EDBBA2D0DC3E94AuL)),
                ((-1, 7, 0x9C7DF410846EC8A1uL, 0xFFA8EC074D09E304uL), (+1, 8, 0xC4277DB81456DC70uL, 0xB86BD9FB64FF03C9uL)),
                ((+1, 9, 0xE9E51B8646E8B1D9uL, 0x2A594B57698D42AAuL), (-1, 7, 0x803EA7630EC3B2E6uL, 0xB0EE4E67057CFA92uL)),
                ((+1, 6, 0xA9199FA27AD4CC9FuL, 0xF513A18EBA71B22CuL), (-1, 7, 0xB22217CDAAFE0B5CuL, 0xFFC1F90F4EE4870AuL)),
                ((-1, 7, 0xF25711B249443AF6uL, 0x934090F2601A0442uL), (+1, 5, 0x88A77144AF62D99FuL, 0x6B4F6F8637DF468AuL)),
                ((+1, 2, 0xB9A5B8D2EA1AD0CFuL, 0xEE6D27A8EF157586uL), (+1, 3, 0x9FDEF1444773C62EuL, 0x80BD031818DFFC79uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p375_0p4375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 0, 0x962D686E6E14D178uL, 0xB63BB842F4E6A0D7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -2, 0x868176133A428698uL, 0x048CAFD58100CEC9uL), (+1, 1, 0xD1C32394415734EEuL, 0x3F27949720162980uL)),
                ((+1, 4, 0x8B347B7DB77C2954uL, 0x69F43A5FA7F2E812uL), (-1, 2, 0x9BA5B6C1CA085348uL, 0x20050C1739E51110uL)),
                ((+1, 3, 0x9E101B56392AE6FAuL, 0x9363A9DBF039F38DuL), (-1, 4, 0x99507E7599710FB0uL, 0xC0379C2C2B53645CuL)),
                ((-1, 6, 0x821E70629C901A59uL, 0x46CE2B74B898020CuL), (+1, 3, 0x8B803EB6783EDFF1uL, 0x39F5418431EDAD73uL)),
                ((-1, 4, 0xDD0D7C1BB2897949uL, 0x7B3FC17DD9B9B5D0uL), (+1, 5, 0x8DDDF5785296C728uL, 0x9A722C01209D1CDFuL)),
                ((+1, 6, 0xAA749EBC5560F442uL, 0x80ED0327C9A883F7uL), (-1, 3, 0x89C61B34A5436723uL, 0xFF82B08E448F9D1EuL)),
                ((+1, 3, 0xF3D827CD0F6DE9D3uL, 0xF3E74537B8B3E66CuL), (-1, 4, 0xA3D6AF6C369AC652uL, 0x3F973DC2D6D2C0BCuL)),
                ((-1, 4, 0xED38F85A55D2B22EuL, 0xBC5CFDD779DD8CFBuL), (+1, 1, 0xDFC3F1866989B234uL, 0xCD3B8FF3CD3C909DuL)),
                (Zero, (+1, 0, 0xC782EE09EDB843A2uL, 0x3354247768483EEDuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p4375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -1, 0xF2C6AF74E34838A2uL, 0x51A847ED1403D688uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0xAA73928382F1F2D7uL, 0x5C4E583060F339E3uL), (+1, 0, 0x82EEFFCB438F06ABuL, 0x46C3B828D788873BuL)),
                ((+1, 3, 0xC62AC8449E417F5EuL, 0xBBBA58DDDA97527FuL), (-1, 3, 0x83D6264379B1A569uL, 0x9B003C101FAF1916uL)),
                ((-1, 4, 0x95D1B9B761CA9F17uL, 0x43F3BDD1DFECD65BuL), (-1, 2, 0xBA1F20973BA38E8AuL, 0x3A2DAB3E3441D42CuL)),
                ((-1, 5, 0xAF7BB0B32D4D2949uL, 0x994EE06A6A93C2B8uL), (+1, 4, 0xB90DA36BEB9504C3uL, 0x580090BF27D01376uL)),
                ((+1, 5, 0xB0363D04DA2367D7uL, 0x3C248EDDC2FB9646uL), (+1, 3, 0x92C6FB941111F947uL, 0x75858E41D5A25C69uL)),
                ((+1, 5, 0xCAE026A15D770E83uL, 0xEAA4B3609B829353uL), (-1, 4, 0xC1BA40FD9C0B39D1uL, 0x702C857BF0F4F399uL)),
                ((-1, 5, 0x8AF6F3B530581C8DuL, 0xBA1F0D1B263BEBD4uL), (-1, 1, 0xDDD0D5179669D8E7uL, 0xBF7FCEC16CDF9375uL)),
                ((-1, 3, 0xDC8885168EA43694uL, 0xB3E6AC26301915E7uL), (+1, 2, 0xD52B4A1EF2732DC0uL, 0xA9F49952928BBD48uL)),
                ((+1, 2, 0x8F5937BCF23458C8uL, 0x4BA563E1BB7F9D32uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0x8C027C97F1540D45uL, 0x0EA7EC7DE67EBBA5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 1, 0xF1492577BA8FF8A1uL, 0x1BB8267F4B14FCE2uL), (+1, 0, 0xC14C350745D37993uL, 0xEA38F5D368A7BC97uL)),
                ((-1, 1, 0xAE83CE98815DD7ABuL, 0xF075A0D5DEE207C3uL), (+1, -1, 0xF49CE780FCD77A90uL, 0xC0BAD675D0BBA6B9uL)),
                ((-1, 0, 0x89EC8470041F6054uL, 0x2F1EB81823D82F15uL), (+1, -2, 0xA859639846B158D9uL, 0x2289D3B23A021562uL)),
                ((-1, -2, 0x81D6CFC5CBE106D5uL, 0x1E56D35AE1E6EE69uL), (+1, -4, 0x88DAB506AE0AA41AuL, 0x547E77ECE9D86B1BuL)),
                ((-1, -5, 0x94C66A18C5955EF4uL, 0x6404FBD01B014EEEuL), (+1, -7, 0x858E2C340DAE7579uL, 0x7A14A348909BAE03uL)),
                ((-1, -9, 0xCB5AE223ED2A5AD1uL, 0xA1C06836D58F3D21uL), (+1, -11, 0x98026A9C04DD2986uL, 0x7491E23D425C528AuL)),
                ((-1, -13, 0x9B70E7B8A9812CB2uL, 0x9F20BA499B0E8CE4uL), (+1, -16, 0xBA2CF4EC31282929uL, 0x81CB74D126B43A39uL)),
                ((-1, -19, 0xE81F75F70840E32FuL, 0xA36739C13C2731EFuL), (+1, -22, 0xCEB382ABB0BB860AuL, 0x84DE213D307CB162uL)),
                ((-1, -26, 0xFA3958D37657B753uL, 0x03225C7E938D01E2uL), (+1, -29, 0x88A296085EB9217DuL, 0x8BE68690794B9311uL)),
                ((-1, -35, 0x81904804FA7CE6CFuL, 0xCEE017EE5E6ECD57uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xA6494D4A96F241A0uL, 0xB743C6E60EBFAAD1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 1, 0xF78B9F979C1C34BCuL, 0xEBA80F9CD79306A8uL), (+1, 0, 0xAC99445B82158B2EuL, 0x1564A6A96DA77A3DuL)),
                ((-1, 1, 0x9FF5C284A340C529uL, 0x96225ECE0C1A3455uL), (+1, -1, 0xC9A25C8A61D3D96CuL, 0xA6ACC9E736B185CBuL)),
                ((-1, -1, 0xEBA4F31B97D1577DuL, 0x67E8496F8909A70CuL), (+1, -2, 0x85CDBDE9833AF71CuL, 0xF7B0A118524EC5D3uL)),
                ((-1, -3, 0xDA6D2EA1841CC2D6uL, 0x79D96555948E34F8uL), (+1, -5, 0xDE607F543DC5DB5CuL, 0xAFF1F949353F36ECuL)),
                ((-1, -5, 0x84BB973C8505993BuL, 0x1270B5420981D3E8uL), (+1, -8, 0xF0991120424D4D6CuL, 0xF6D055753785E007uL)),
                ((-1, -9, 0xD5D13267464FCE02uL, 0x2F8CCB303E9AF722uL), (+1, -11, 0xAAD2BDE307252C7EuL, 0x145C3440FE3024C9uL)),
                ((-1, -13, 0xE1B94E106CEBDA0AuL, 0x5FFA4E753754B6D7uL), (+1, -15, 0x9CC13902F2C9B449uL, 0xCA9395C102B4C41EuL)),
                ((-1, -17, 0x971583C34F3DE552uL, 0x4458079FE572BDCDuL), (+1, -20, 0xB2B1F7E8F4465CE6uL, 0x88BE3711DD3392D3uL)),
                ((-1, -23, 0xF12A0C13E1A167F8uL, 0x20D8E469C3D8F92EuL), (+1, -26, 0xEB485DCC63ECFE3DuL, 0x001C612C81812F9AuL)),
                ((-1, -29, 0xCE06F86EB90C9E10uL, 0x0E936449DA364BB4uL), (+1, -32, 0x9D01A6E5C3328A52uL, 0x6F5D9E62AC420E39uL)),
                ((-1, -36, 0x99AF1D8A0AB4773AuL, 0x09C391D3A0E5C86AuL), (+1, -40, 0xA3A09F66EDA307ADuL, 0xEC4D7B3D123D3B27uL)),
                ((-1, -46, 0xF69D65ED08688BD2uL, 0x4B17BC8447B6CB5DuL), (+1, -50, 0x80A60235EAAF06EDuL, 0xE7AEBEC93E5C9C09uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xEB1B56D761E029F9uL, 0x532E118565F0111DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 1, 0xAB32CC7E32CA6C4FuL, 0x5CB59960B74125BAuL), (+1, -1, 0xABF336DD7187BB01uL, 0x47026D39B244F01DuL)),
                ((-1, -1, 0xD8D505486A831F11uL, 0x241163B4BDC4ABD5uL), (+1, -3, 0xC8146ABC523250B9uL, 0xE93237D11497ABCDuL)),
                ((-1, -3, 0x9CC4671E1669ECC8uL, 0x54AD324067BCADCAuL), (+1, -5, 0x8437031F54C14918uL, 0x1A53B50FD2D713D9uL)),
                ((-1, -6, 0x8ED4B7C2AF7A2D90uL, 0x29A8BD50950988CAuL), (+1, -9, 0xDAC212BFFCB67EDEuL, 0x812A03D96EDA0EB1uL)),
                ((-1, -10, 0xAAD689D6CE67FCF9uL, 0x159D959876A8D309uL), (+1, -13, 0xEB8FECCBDB63B367uL, 0xCA153858D9533664uL)),
                ((-1, -14, 0x87925D511BBE91DCuL, 0x50EA156E29F3E7C9uL), (+1, -17, 0xA66795D948A26D29uL, 0xFE86F69990C3F3D4uL)),
                ((-1, -19, 0x8D2877BF8C6DB0BFuL, 0x66147329D898F0EDuL), (+1, -22, 0x97E0F9BDA345BCE6uL, 0x5CC8E0F061B4A8E3uL)),
                ((-1, -25, 0xBA8D73894FE4EFB6uL, 0x46ACF5081CE99C21uL), (+1, -28, 0xAC246B550BF9339CuL, 0xAA9F9636BD879B31uL)),
                ((-1, -31, 0x93209B51CAC2A53CuL, 0x9161C7705BE1AA7FuL), (+1, -35, 0xE1434B5CFEE27F77uL, 0x95C77439D824CC8FuL)),
                ((-1, -39, 0xF8A20E57F1767E7BuL, 0xA43E2F6F0452C85DuL), (+1, -42, 0x9554F9D9AFC453D8uL, 0x7E4943988DB31BA8uL)),
                ((-1, -47, 0xB79B53AAA028A4E1uL, 0x4F2E7F3CEB2F9F02uL), (+1, -51, 0x9A87A76BF748C9C2uL, 0x4B3275F6850D3055uL)),
                ((-1, -57, 0x91FC57A769C22ECDuL, 0x3861D390BA2D9B0EuL), (+1, -63, 0xF11EA2AC1C6F378BuL, 0x9EB4CDE79C232E3CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 2, 0x9D8006B35C13472EuL, 0xC68BBF1ADC1B54A2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xE2582E2C6F61B2C4uL, 0x226413BA7FEDFDD8uL), (+1, -2, 0xAB3706B125D66936uL, 0xFC19AECAE6EA1134uL)),
                ((-1, -2, 0x8D8900B060EF6359uL, 0x9857E9D739781034uL), (+1, -5, 0xC657E480269B694AuL, 0xAD6C78622D9FA4D4uL)),
                ((-1, -6, 0xCA33D19A6DA46292uL, 0x7B216C46E0F7A51DuL), (+1, -8, 0x8277A26B19B1D439uL, 0xAEDB0EC687A7E86CuL)),
                ((-1, -10, 0xB61EFF88D1A71416uL, 0x42F56D59518FEFFCuL), (+1, -13, 0xD6D7B7A4C812F505uL, 0xF51BEEAD8C989DBDuL)),
                ((-1, -15, 0xD775A6EAF254A3F8uL, 0x4312978CD557D19BuL), (+1, -18, 0xE63568CEC9C5D14BuL, 0x2D69C79C0A935700uL)),
                ((-1, -20, 0xA9357DE452F928D9uL, 0x81F5DA8A82B07050uL), (+1, -23, 0xA1CB4DBD07015CD6uL, 0xAD8B3DD6197F2EC7uL)),
                ((-1, -26, 0xAE71FEC6B4887DA7uL, 0x11BDC51CEFAE52F5uL), (+1, -29, 0x92E424552A3DB7CAuL, 0x0C458C34A4627C07uL)),
                ((-1, -33, 0xE464551A5D29E6CAuL, 0xE6454D77E28CA2B4uL), (+1, -36, 0xA59440535444C9F2uL, 0x8D5B888E9E5C9435uL)),
                ((-1, -40, 0xB289A9B2445E6066uL, 0x69F4C17CF2F6154AuL), (+1, -44, 0xD773013231E77471uL, 0xC6528F1ADF01FAA9uL)),
                ((-1, -48, 0x959C4E74B88F89A4uL, 0xAE851319B75E4064uL), (+1, -52, 0x8DFDD7D3400B27C7uL, 0xC99040AEDD0185FDuL)),
                ((-1, -58, 0xDB4473FD8D60D37BuL, 0xE1E086B18138D425uL), (+1, -62, 0x920BD9E2E0D6457DuL, 0x2F777526A4E6D554uL)),
                ((-1, -69, 0xAD1D078ABBA5216DuL, 0x4C819F77838C62ACuL), (+1, -75, 0xE272F44FBC9B3E6BuL, 0x4328B0C556582BF0uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 2, 0xCD430E4046E7DC62uL, 0xA8DDF51920A19CBBuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0x923A02F863F866FBuL, 0x1177A20606163F65uL), (+1, -3, 0xAA8F11219D4EB9E8uL, 0xFC820770EE4EDBB0uL)),
                ((-1, -4, 0xB55AE31FD3EF48A0uL, 0x24B206E4EA163138uL), (+1, -7, 0xC4CEB8D7AC719D01uL, 0xD13B61067587F77EuL)),
                ((-1, -8, 0x80812DBCD0D13C5BuL, 0x8DD7C8E0CBAAFDA7uL), (+1, -11, 0x80EFB4A42BD1ACA9uL, 0x656A7179D8C6D624uL)),
                ((-1, -14, 0xE5B013FE9E6AC953uL, 0xB2C289476EB8CED3uL), (+1, -17, 0xD3726C6BA9034FD7uL, 0xA0BAC9579EE1F0B4uL)),
                ((-1, -19, 0x86D963C6743165A6uL, 0x413C8187D3D41D18uL), (+1, -23, 0xE19D6FD8EF913BB5uL, 0xF1F9B442FE4D6717uL)),
                ((-1, -26, 0xD2470580E5488A16uL, 0xEBAC393FB124835EuL), (+1, -29, 0x9DE1E15123C875A6uL, 0x62C90A1276E2E655uL)),
                ((-1, -33, 0xD748B021670B7E86uL, 0x025B1C26AAD5351FuL), (+1, -36, 0x8EB58892FF99A598uL, 0x694521688F7EA114uL)),
                ((-1, -40, 0x8BFF335E1EC91BF9uL, 0xD56A2ED74592F302uL), (+1, -44, 0xA024D893C64DDC4DuL, 0x8F4E2A572C744374uL)),
                ((-1, -49, 0xD97F31945E9209A0uL, 0x3EA9379ADF33FC13uL), (+1, -53, 0xCF6D70AF0722EAFCuL, 0x1E674DCB6064FE87uL)),
                ((-1, -58, 0xB52C000F981EA643uL, 0x57807243F9CEE638uL), (+1, -62, 0x8812B563A66DA45CuL, 0x0C21A87C86F26D00uL)),
                ((-1, -68, 0x8405424762A2E48DuL, 0xDB001F61999976B3uL), (+1, -73, 0x8B4DD69CFB8BB97CuL, 0x4E9337F14FE4BD97uL)),
                ((-1, -81, 0xCF62D2C77B769E85uL, 0x97F82B403662B7BBuL), (+1, -87, 0xD6F8474C4B9D0294uL, 0xA1D356CEE8B62B52uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 3, 0x83C297D0B98E0291uL, 0x97541B23E165C163uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -1, 0xC21281E75E04C21CuL, 0xD7777C46A0737554uL), (+1, -4, 0xB13A9B864B0BB0E0uL, 0xB0958404A8F2D0F1uL)),
                ((-1, -6, 0xFA469EA652B806E4uL, 0xB90A1FEF882BC13DuL), (+1, -9, 0xD5CA5897C23411F4uL, 0x9DF9508F598B01D2uL)),
                ((-1, -11, 0xB9B6EE003034FD49uL, 0x2879BC3130D39B23uL), (+1, -14, 0x938580EDA8CDF2F3uL, 0xD47E4A3F26C80FD5uL)),
                ((-1, -17, 0xAF5A556D7968FFF0uL, 0xF1337742B4D86061uL), (+1, -20, 0x809B9366CD6AEE04uL, 0xE00F4523120BBE88uL)),
                ((-1, -24, 0xDBFD646B9DF6D1F4uL, 0xD2D92465768345E2uL), (+1, -27, 0x93AAE0152588C6ABuL, 0x6F6611FEEC25F525uL)),
                ((-1, -31, 0xB9F326720AEA43A5uL, 0xA07A69193CDF7B2AuL), (+1, -35, 0xE1F3BF5FCEA1CA88uL, 0xF9C09E5D48F88096uL)),
                ((-1, -39, 0xD2757A8C5A5F3C43uL, 0x3A37F4CA208A87CAuL), (+1, -43, 0xE426467FD784347DuL, 0x8FF6A9DDC89BC7BEuL)),
                ((-1, -47, 0x9B76B905C2B5428DuL, 0x0628954EF0A6733BuL), (+1, -51, 0x936EAC0191E7FA04uL, 0x1A634A56E3FCA990uL)),
                ((-1, -56, 0x8EC718104C134BD8uL, 0x53291532A4375E1BuL), (+1, -61, 0xE64E6CA5011A016DuL, 0x2E47ABF01B048084uL)),
                ((-1, -66, 0x95E66829EBC284CDuL, 0x6FBCD5A56A8FF41BuL), (+1, -71, 0xC4B7B5AC25B6D820uL, 0x898785426286B300uL)),
                ((-1, -77, 0x9AD85EC6607CFC17uL, 0x2DE8D9BBCC1BF3B4uL), (+1, -82, 0x9882B4C8803CBCEFuL, 0xA90FBBA7D2C4093CuL)),
                ((-1, -90, 0xE5E94B97698C4A73uL, 0xA00E048EC5129258uL), (+1, -95, 0x8B3ED4E2E9A9CAAAuL, 0x4DD2E2F6BD36642DuL)),
                ((-1, -105, 0xA1C4C662E93FF61CuL, 0xAC962DBFD36EFDD4uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm128_256 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 3, 0xA7C0BEE83773C451uL, 0x36C6E754D498E570uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -2, 0xF62D00F1C5479DE6uL, 0x797BC7A7573224CAuL), (+1, -5, 0xB0D13C431DFA0B59uL, 0x22A4B5AD2B10B9FCuL)),
                ((-1, -7, 0x9E2907AA09DB6144uL, 0x8CF33A7A6D285F50uL), (+1, -11, 0xD4CB0C1A364DCB25uL, 0x2561D043AE4F6374uL)),
                ((-1, -14, 0xE9E4184BBDDA1DCFuL, 0xF1631D6C053B5ED7uL), (+1, -17, 0x927C386B057B4BF1uL, 0x0692A57F93DA3D33uL)),
                ((-1, -21, 0xDC139A74C4C599ADuL, 0x10E7373767B54E45uL), (+1, -25, 0xFECC3C8B3D6A7E8DuL, 0x365048F2B387F14DuL)),
                ((-1, -28, 0x899525D35AB34E2FuL, 0xAA49E067DBDF6A1EuL), (+1, -32, 0x91ED57F5B2344F4DuL, 0x58A4A14ABDF71578uL)),
                ((-1, -37, 0xE7D29C37F69209B0uL, 0xC5B8EEA541F1D5A4uL), (+1, -41, 0xDEBF8DA8C3673B84uL, 0x4A0B20687201ECD3uL)),
                ((-1, -45, 0x82C5332272828EE0uL, 0x759F3C01797973B4uL), (+1, -50, 0xE05E279E11DB0A5CuL, 0x4F69A227A51DA1DDuL)),
                ((-1, -55, 0xC099DE4AFC58959FuL, 0x538A2B3544CDF288uL), (+1, -59, 0x90A2C5339B3A60FCuL, 0x9FDA8DDB4351C541uL)),
                ((-1, -65, 0xB05BE85636B9123AuL, 0x7B9622B5AB1DE9ADuL), (+1, -70, 0xE163B4513F05FE2DuL, 0x2BBBB789C28981DDuL)),
                ((-1, -76, 0xB8A134EDA0A01E32uL, 0xC6352D236A60AE7EuL), (+1, -81, 0xC00D9BD99464FDE4uL, 0xD750977E65A7986CuL)),
                ((-1, -88, 0xBE3450DFF855695DuL, 0x145BC12D9C1B9114uL), (+1, -93, 0x9489CF24AAA452C2uL, 0xF63E116AA1FC7B38uL)),
                ((-1, -101, 0x8CD864B002EF29D5uL, 0xC57A2B19CFABA6EEuL), (+1, -107, 0x874C702A76A05E3AuL, 0xDF35B2A52FA53D31uL)),
                ((-1, -118, 0xC5BAC4ADB7EFBB73uL, 0x240DB1866B2BC885uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm256_512 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 3, 0xD495D7C3D3DB8D6FuL, 0xFBE62894AD6390B1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -2, 0x9B987834C94CDD38uL, 0xE8B170E09E0F6C99uL), (+1, -6, 0xB0823E07CF9C0D0AuL, 0xCB128AEB2E91C8D1uL)),
                ((-1, -9, 0xC771224199E75F25uL, 0x7A93D3E3EE70238EuL), (+1, -13, 0xD40C8A5CCCB9E8EFuL, 0xC6C4A78C383BBE7EuL)),
                ((-1, -16, 0x931DA62AD2E0AACAuL, 0x8DE68B5AB32D3353uL), (+1, -20, 0x91B72FCF53ACBA91uL, 0x1C653519A1B095F0uL)),
                ((-1, -24, 0x8A19EA61DE4386BCuL, 0x201BE0CCA8A5D420uL), (+1, -29, 0xFD02B98CD2D9BD18uL, 0xE52021342AD34ED3uL)),
                ((-1, -33, 0xAC45F407FA7C5631uL, 0x8B25FB35E2B197EFuL), (+1, -37, 0x90A593DF04E7A1DCuL, 0x1F51CB70C3D5CB29uL)),
                ((-1, -42, 0x90CFF3E405260CBBuL, 0x47841BBB611205D2uL), (+1, -47, 0xDC67153BAC358CB8uL, 0x050EA9DE90FE0715uL)),
                ((-1, -52, 0xA3051B663374316FuL, 0x38A41649B41C1B33uL), (+1, -57, 0xDD9CDB7718EDC8A3uL, 0x19193592FA8A31C4uL)),
                ((-1, -63, 0xEF97759A981E9090uL, 0xA2AC4F219E27188DuL), (+1, -67, 0x8E9BCBEC3480411BuL, 0xF6E9B160317A69E9uL)),
                ((-1, -74, 0xDAF0893E1AF3135AuL, 0x004D3A2307E23848uL), (+1, -79, 0xDDD7C4D674124A3FuL, 0x188648197B618FEBuL)),
                ((-1, -86, 0xE4C22F5ADF46C808uL, 0xA0C90B65ABE989ECuL), (+1, -91, 0xBCB49836B0641978uL, 0x8238C4AF6A1918D4uL)),
                ((-1, -99, 0xEB39B44DD78063F1uL, 0x5D19BE3D34F048F8uL), (+1, -104, 0x91B3A5C9D44FC8EDuL, 0x380B934C423149DFuL)),
                ((-1, -113, 0xADDFFFFAAF619755uL, 0x2E3D8C98D95EE0BEuL), (+1, -119, 0x847E951195AFC121uL, 0x000A5198A14E0252uL)),
                ((-1, -131, 0xF3B1DC26FE94727AuL, 0x8BEEF8F150CD9958uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm512_1024 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 4, 0x8659892AFFD04B69uL, 0x851ADA726D80FA6DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -3, 0xC4561CF4795A8A50uL, 0x908FE619FE96D2F1uL), (+1, -7, 0xB0484A52DC445FF3uL, 0x46EAA29993D7C8E9uL)),
                ((-1, -11, 0xFB3F4CFFC11F0732uL, 0x58E67F454615929EuL), (+1, -15, 0xD3814AF675465468uL, 0x350800C83BCD5F36uL)),
                ((-1, -19, 0xB90737B74C4F05FEuL, 0x771B3EB16F4FE467uL), (+1, -23, 0x9127B3C001CDFD7FuL, 0x56BDBE11DFAC808EuL)),
                ((-1, -28, 0xAD69DE3842627A4FuL, 0x19EDEC1214F58268uL), (+1, -33, 0xFBB6CBBE914269E1uL, 0x4B4A3C922DECA4F7uL)),
                ((-1, -38, 0xD7FC4B60B5AEFAA2uL, 0xD075A4C4CEA9344EuL), (+1, -42, 0x8FB8ABDFAF346914uL, 0x6533812ED5CAD0A6uL)),
                ((-1, -48, 0xB5479E11C10D42A8uL, 0x30D779A87AA436C0uL), (+1, -53, 0xDAB6B373D42DABF1uL, 0x31916A4796880477uL)),
                ((-1, -59, 0xCBC4F1DA1A272BE7uL, 0x48CA6B6B335F79E7uL), (+1, -64, 0xDBA2E858D21B7DB1uL, 0x26F65E759D2051A9uL)),
                ((-1, -70, 0x9586501CAACB7690uL, 0x3759A6D456453964uL), (+1, -75, 0x8D28E9B49D5AC8B3uL, 0x2FDB2E88BC7BE86EuL)),
                ((-1, -82, 0x887210DD291135CBuL, 0x469B7C8CA0B1AEA6uL), (+1, -88, 0xDB515CD3D3B8CECFuL, 0xDC685815B1BFBDB0uL)),
                ((-1, -95, 0x8E5FA28ACC6BEA32uL, 0xC29928D3EF441DDFuL), (+1, -101, 0xBA54BD372BE2B139uL, 0x7F74340B0224C95CuL)),
                ((-1, -109, 0x9235E47A35BA5EF2uL, 0xDA3F19A8F70F14C6uL), (+1, -115, 0x8FB2893CDAEA6E2FuL, 0x56BB03A600C42946uL)),
                ((-1, -125, 0xD7E3BF7A04D84F96uL, 0x5A6541FECCE7291AuL), (+1, -131, 0x82853E5899BBD4D4uL, 0xF1CFA7DA937130D8uL)),
                ((-1, -143, 0x971DE2C241D784ECuL, 0xE506B71AABB9256DuL), Zero),
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
                    else {
                        y = ApproxUtil.Pade(x - 0.375d, pade_upper_0p375_0p5);
                    }

                    return y;
                }
                else {
                    ddouble v;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 3)), pade_upper_expm3_4);
                    }
                    else if (exponent >= -6) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 4)), pade_upper_expm4_6);
                    }
                    else if (exponent >= -8) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 6)), pade_upper_expm6_8);
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
                    else if (exponent >= -80) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 64)), pade_upper_expm64_80);
                    }
                    else {
                        v = 1d / Cbrt(Ldexp(PI, 1));
                    }

                    ddouble y = v / ExMath.Pow2d3(x);

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
                    else if (x <= 0.4375d) {
                        y = ApproxUtil.Pade(x - 0.375d, pade_lower_0p375_0p4375);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.4375d, pade_lower_0p4375_0p5);
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
