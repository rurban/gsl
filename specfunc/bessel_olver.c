/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_airy.h"
#include "gsl_sf_chebyshev.h"

#define CubeRoot2_ 1.25992104989487316476721060728


/* Chebyshev fit for f(x) = z(x)^6 A_3(z(x)),  z(x) = 22/(10(x+1)) */
static double A3_gt1_data[31] = {
  -0.123783199829515294670493131190,
   0.104636462534700704670877382304,
  -0.067500816575851826744877535903,
   0.035563362418888483652711005520,
  -0.0160738524035979408472979609051,
   0.0064497878252851092073278056238,
  -0.00235408261133449663958121821593,
   0.00079545702851302155411892534965,
  -0.000252149207458550798957848256370,
   0.000075740045960693929211533018331,
  -0.0000217291796633962343440797826255,
   5.9914810727868915476543145465e-6,
  -1.59587815718089921629537198165e-6,
   4.1232986512903717525448312012e-7,
  -1.03697259934176591019139191006e-7,
   2.54579823042665411459992350218e-8,
  -6.1161715053791743082427422443e-9,
   1.44093461991386588878714613196e-9,
  -3.3350445956255561668232014995e-10,
   7.5950686572918996453336138108e-11,
  -1.70422963344094303773899002780e-11,
   3.7723525020626230919721640081e-12,
  -8.2460237635733980528416501227e-13,
   1.78169615279977976962518688754e-13,
  -3.8084101506541792942694560802e-14,
   8.0593669930916099079755351563e-15,
  -1.68965659616417390174526369644e-15,
   3.5115651805888443184822853595e-16,
  -7.2384771938569255638904297651e-17,
   1.48065989776771761062838402443e-17,
  -3.00692857507873036348979979631e-18
};
static struct gsl_sf_cheb_series A3_gt1_cs = {
  A3_gt1_data,
  30,
  -1,1,
  (double *)0,
  (double *)0
};

/* chebyshev expansion for f(x) = z(x)^8 A_4(z(x)), z(x) = 12/(5(x+1)) */
static double A4_gt1_data[30] = {
  1.15309329391198493586724229008,
 -1.01812701728669338904729927846,
  0.71964022270555684403652781941,
 -0.42359963977172689685150061355,
  0.215024488759339557817435404261,
 -0.096751915348145944032096342479,
  0.039413982058824310099856035361,
 -0.0147752256925616979637811150142,
  0.0051621145141593705169478232709,
 -0.00169783446445524322560925166335,
  0.00052995667873006847211519193478,
 -0.000158020275749964771156679748555,
  0.000045254366680989687988902825193,
 -0.0000125037229654746380154886009674,
  3.3457656998119148699124716204e-6,
 -8.6981575241150758412492331833e-7,
  2.20308954843256456408239406250e-7,
 -5.4493369492600677068285936533e-8,
  1.31904572817248291071393855558e-8,
 -3.13015601833773791589511917685e-9,
  7.2937802527123344842593076131e-10,
 -1.67120801379451404073489401094e-10,
  3.7700053248213600430503521194e-11,
 -8.3824538848817227637828899571e-12,
  1.83887419100497668652740371944e-12,
 -3.9835919980753778560117573063e-13,
  8.5288827136546615604290389711e-14,
 -1.80602278691144169986532668360e-14,
  3.7849342199690728470461022877e-15,
 -7.8552867468122209577151823365e-16
};
static struct gsl_sf_cheb_series A4_gt1_cs = {
  A4_gt1_data,
  17, /* 29, */
  -1, 1,
  (double *)0,
  (double *)0
};

/* Chebyshev fit for f(x) = z(x)^3 B_2(z(x)), z(x) = 12/(5(x+1)) */
static double B2_gt1_data[40] = {
  0.00118587147272683864479328868589,
  0.00034820459990648274622193981840,
 -0.000304113044256397681030758645670,
  0.0000281206628401234353148468288586,
  0.000044935252959016131844898987478,
 -0.0000303762999709307219677948967641,
  0.0000112597964712387572194974397039,
 -2.48325339695177759919510082181e-6,
 -9.9003813640537799587086928278e-8,
  4.9259859656183110299492296029e-7,
 -3.7644120964426705960749504975e-7,
  2.28878285213346251896391225087e-7,
 -1.32026873708222037314898550504e-7,
  7.7019669092537400811434860763e-8,
 -4.6589706973010511603890144294e-8,
  2.93964762330139237119785229628e-8,
 -1.92932306119882829191019545384e-8,
  1.30991070137287178424069068964e-8,
 -9.1509111940885962831104149355e-9,
  6.5483472971925614347299375295e-9,
 -4.7831253582139967461241674569e-9,
  3.5562625457426178152760148639e-9,
 -2.68533894440084141869165621038e-9,
  2.05547386671342001457818572890e-9,
 -1.59231720195174262778865227580e-9,
  1.24659232134643814573194814977e-9,
 -9.8494846881180588507969988989e-10,
  7.8438674499372126663957464312e-10,
 -6.2877567918342950225937136855e-10,
  5.0662318868755257959686944117e-10,
 -4.0962270881243451160378710952e-10,
  3.3168684677374908553161911299e-10,
 -2.68294066198474506335961633057e-10,
  2.16039881221845683755610778737e-10,
 -1.72323733095602784020121244819e-10,
  1.35127090896114706266178304344e-10,
 -1.02853547325386630131675797924e-10,
  7.4211345443901713467637018423e-11,
 -4.8124980266864320351456993068e-11,
  2.36665346944763060774168319576e-11
};
static struct gsl_sf_cheb_series B2_gt1_cs = {
  B2_gt1_data,
  39,
  -1, 1,
  (double *)0,
  (double *)0
};


/* Chebyshev fit for f(x) = z(x)^6 B_3(z(x)), z(x) = 12/(5(x+1)) */
static double B3_gt1_data[30] = {
 -0.0102445379362695740863663926486,
  0.0036618484329295342954730801917,
  0.00261542524985993032825693211170,
 -0.0036187389410353156728771706336,
  0.00218785641576922759446134524620,
 -0.00082199523035908035844265168213,
  0.000128177388915563149432131652022,
  0.000100094465336803298572054863717,
 -0.000128829334466377427345314778787,
  0.000101362642026965138678214872049,
 -0.000070002758496595562219165727328,
  0.000046948863967574304316079551463,
 -0.000031900038697178376863569456964,
  0.0000223145366844777521966594747879,
 -0.0000161110219771243953930033643785,
  0.0000119663442499073521446663351298,
 -9.0986920398931223804111374679e-6,
  7.0492613694235423068926562567e-6,
 -5.5425216624642184684300615394e-6,
  4.4071884714230296614449244106e-6,
 -3.5328595506791663127928952625e-6,
  2.84594975572077091520522824686e-6,
 -2.29592697828824392391071619788e-6,
  1.84714740375289956396370322228e-6,
 -1.47383331248116454652025598620e-6,
  1.15687781098593231076084710267e-6,
 -8.8174688524627071175315084910e-7,
  6.3705856964426840441434605593e-7,
 -4.1358791499961929237755474814e-7,
  2.03541511587388198674779968076e-7
};
static struct gsl_sf_cheb_series B3_gt1_cs = {
  B3_gt1_data,
  29,
  -1, 1,
  (double *)0,
  (double *)0
};


/* Chebyshev fit for f(x) = z(x) B_2(z(x)), z(x) = 2(x+1)/5 */
static double B2_lt1_data[40] = {
  0.00073681565841337130021924199490,
  0.00033803599647571227535304316937,
 -0.000082517232192397540242105526788,
 -0.000033908799486564325459007797102,
  0.0000196139805684888181669401488816,
 -2.35593745904151401624656805567e-6,
 -1.79055017080406086541563835433e-6,
  1.33129571185610681090725934031e-6,
 -5.3887944471543654413067395617e-7,
  1.49603056041381416881299945557e-7,
 -1.83377228267274327911131293091e-8,
 -1.33191430762944336526965187651e-8,
  1.60642096463700438411396889489e-8,
 -1.28932576330421806740136816643e-8,
  9.6169275086179165484403221944e-9,
 -7.1818502280703532276832887290e-9,
  5.4744009217215145730697754561e-9,
 -4.2680446690508456935030086136e-9,
  3.3941665009266174865683284781e-9,
 -2.74407140722216738821631351700e-9,
  2.24883615221082552291930389629e-9,
 -1.86382407166087488620879233377e-9,
  1.55923509408053735008664404016e-9,
 -1.31457439377323306092426330705e-9,
  1.11537167772150478427902449682e-9,
 -9.5117576805266622854647303110e-10,
  8.1428799553234876296804561100e-10,
 -6.9893770813548773664326279169e-10,
  6.0073113636087448745018831981e-10,
 -5.1627434258513453901420776514e-10,
  4.4290993195074905891788459756e-10,
 -3.7852978599966867611179315200e-10,
  3.2143959338863177145307610452e-10,
 -2.70259266806207775949922211438e-10,
  2.23848577724579185392282343216e-10,
 -1.81250716642766780465512717011e-10,
  1.41648700087136687672930085466e-10,
 -1.04331018571327824858133259817e-10,
  6.8663910168392483929411418190e-11,
 -3.4068313177952244040559740439e-11
};
static struct gsl_sf_cheb_series B2_lt1_cs = {
  B2_lt1_data,
  39,
  -1, 1,
  (double *)0,
  (double *)0
};


/* Chebyshev fit for f(x) = B_3(2(x+1)/5) */
static double B3_lt1_data[40] = {
 -0.00137160820526992057354001614451,
 -0.000254749379511010499826805613016,
  0.000247629755478958816520734677711,
  0.000052296572814801967493139302652,
 -0.000074883542726215123850165937596,
  0.0000141688001289104644998044974579,
  0.0000152898606017218369074257622919,
 -0.0000166867229707859051429332532645,
  0.0000106176518953645901873958509440,
 -5.8220577442406209989680801335e-6,
  3.3322423743855900506302033234e-6,
 -2.23292405803003860894449897815e-6,
  1.74816651036678291794777245325e-6,
 -1.49581306041395051804547535093e-6,
  1.32759146107893129050610165582e-6,
 -1.19376077392564467408373553343e-6,
  1.07878303863211630544654040875e-6,
 -9.7743335011819134006676476250e-7,
  8.8729318903693324226127054792e-7,
 -8.0671146292125665050876015280e-7,
  7.3432860378667354971042255937e-7,
 -6.6897926072697370325310483359e-7,
  6.0966619703735610352576581485e-7,
 -5.5554095284507959561958605420e-7,
  5.0588335673197236002812826526e-7,
 -4.6008146297767601862670079590e-7,
  4.1761348515688145911438168306e-7,
 -3.7803230006989446874174476515e-7,
  3.4095248501364300041684648230e-7,
 -3.06039597513547495206150154725e-7,
  2.73001341793656905896404589934e-7,
 -2.41580282507623047560442542315e-7,
  2.11547810382987519856891138686e-7,
 -1.82699113287567712014652233136e-7,
  1.54848950858085137490261730740e-7,
 -1.27828068515558093692264404959e-7,
  1.01480117253948925651742073418e-7,
 -7.5658969771439627809239950461e-8,
  5.0226342286491286957075289622e-8,
 -2.50496456602598829705475558313e-8
};
static struct gsl_sf_cheb_series B3_lt1_cs = {
  B3_lt1_data,
  39,
  -1, 1,
  (double *)0,
  (double *)0
};


/* Chebyshev fit for f(x) = A_3(9(x+1)/20) */
static double A3_lt1_data[40] = {
  -0.000179825614721344185876349801172,
  -0.00036558603837525275836608884064,
  -0.0000281939805592962885029440636281,
   0.000167045398638757367698127860671,
  -0.000070989699703476743076230448503,
  -8.4470843942344237748899879940e-6,
   0.0000273413090343147765148014327150,
  -0.0000199073838489821681991178018081,
   0.0000100004176278235088881096950105,
  -3.9739852013143676487867902026e-6,
   1.22653577664495743068826932667e-6,
  -1.88755584306424047416914864854e-7,
  -1.37482206060161206336523452036e-7,
   2.10326379301853336795686477738e-7,
  -2.05583778245412633433934301948e-7,
   1.82377384812654863038691147988e-7,
  -1.58130247846381041027699152436e-7,
   1.36966982725588978654041029615e-7,
  -1.19250280944620257443805710485e-7,
   1.04477169029350256435316644493e-7,
  -9.2064832489437534542041040184e-8,
   8.1523798290458784610230199344e-8,
  -7.2471794980050867512294061891e-8,
   6.4614432955971132569968860233e-8,
  -5.7724095125560946811081322985e-8,
   5.1623107567436835158110947901e-8,
  -4.6171250746798606260216486042e-8,
   4.1256621998650164023254101585e-8,
  -3.6788925543159819135102047082e-8,
   3.2694499457951844422299750661e-8,
  -2.89125899697964696586521743928e-8,
   2.53925288725374047626589488217e-8,
  -2.20915707933726481321465184207e-8,
   1.89732166352720474944407102940e-8,
  -1.60058977893259856012119939554e-8,
   1.31619294542205876946742394494e-8,
  -1.04166651771938038563454275883e-8,
   7.7478015858156185064152078434e-9,
  -5.1347942579352613057675111787e-9,
   2.55835415945867239672615043215e-9
};
static struct gsl_sf_cheb_series A3_lt1_cs = {
  A3_lt1_data,
  39,
  -1, 1,
  (double *)0,
  (double *)0
};

/* chebyshev fit for f(x) = A_4(2(x+1)/5) */
static double A4_lt1_data[30] = {
  0.000090547037700516109469582267356,
  0.00033066000498098017589672988293,
  0.000197374537343639891272260732721,
 -0.000154908097259320377200347628890,
 -0.000045149489355387300854792804540,
  0.000079768817826039408894445739238,
 -0.000033145661545447409862649932513,
 -1.88212148790135672249935711657e-6,
  0.0000114788756505519986352882940648,
 -9.2263039911196207101468331210e-6,
  5.1401128250377780476084336340e-6,
 -2.38418218951722002658891397905e-6,
  1.00664292214481531598338960828e-6,
 -4.2322467809649006026424997054e-7,
  2.00132031535793489976535190025e-7,
 -1.18689501178886741400633921047e-7,
  8.7819524319114212999768013738e-8,
 -7.3964150324206644900787216386e-8,
  6.5780431507637165113885884236e-8,
 -5.9651053193022652369837650411e-8,
  5.4447762662767276209052293773e-8,
 -4.9802057381568863702541294988e-8,
  4.5571368194694340198117635845e-8,
 -4.1682117173547642845382848197e-8,
  3.8084701352766049815367147717e-8,
 -3.4740302885185237434662649907e-8,
  3.16165570647015106112736920601e-8,
 -2.86857394876895562523748792670e-8,
  2.59237521171322544290027965998e-8,
 -2.33094285521905873046628834773e-8
};
static struct gsl_sf_cheb_series A4_lt1_cs = {
  A4_lt1_data,
  29,
  -1, 1,
  (double *)0,
  (double *)0
};


/* checked OK [GJ] Thu Apr 30 22:41:56 MDT 1998 */
static double olver_B0(double z, double abs_zeta)
{
  if(z < 0.98) {
    double t = 1./sqrt((1-z)*(1+z));
    return -5./(48.*abs_zeta*abs_zeta) + t*(-3 + 5.*t*t)/(24.*sqrt(abs_zeta));
  }
  else if(z < 1.02) {
    double a = 1.-z;
    return  0.0179988721413553309252458658183
           +0.0111992982212877614645974276203  * a
	   +0.0059404069786014304317781160605  * a*a
	   +0.0028676724516390040844556450173  * a*a*a
	   +0.0012339189052567271708525111185  * a*a*a*a
	   +0.0004169250674535178764734660248  * a*a*a*a*a
	   +0.00003301733850859498069527773655 * a*a*a*a*a*a
	   -0.00013180762385782030099901064251 * a*a*a*a*a*a*a
	   -0.00019068703700508472398139456474 * a*a*a*a*a*a*a*a
	   ;
  }
  else {
    double t = 1./sqrt((z-1)*(z+1));
    return -5./(48.*abs_zeta*abs_zeta) + t*( 3 + 5.*t*t)/(24.*sqrt(abs_zeta));
  }
}

/* checked OK [GJ] Thu Apr 30 21:50:09 MDT 1998 */
static double olver_B1(double z, double abs_zeta)
{
  if(z < 0.88) {
    double t   = 1./sqrt((1-z)*(1+z));
    double t2  = t*t;
    double rz  = sqrt(abs_zeta);
    double z32 = rz*rz*rz;
    double z92 = z32*z32*z32;
    double term1 = t*t*t * (30375. - 369603.*t2 + 765765.*t2*t2 - 425425.*t2*t2*t2)/414720.;
    double term2 = 85085./(663552.*z92);
    double term3 = 385./110592.*t*(3.-5.*t2)/(abs_zeta*abs_zeta*abs_zeta);
    double term4 = 5./55296.*t2*(81. - 462.*t2 + 385.*t2*t2)/z32;
    return -(term1 + term2 + term3 + term4)/rz;
  }
  else if(z < 1.12) {
    double a = 1.-z;
    return -0.00149282953213429172050073403334
           -0.00175640941909277865678308358128 * a
	   -0.00113346148874174912576929663517 * a*a
	   -0.00034691090981382974689396961817 * a*a*a
	   +0.00022752516104839243675693256916 * a*a*a*a
	   +0.00051764145724244846447294636552 * a*a*a*a*a
	   +0.00058906174858194233998714243010 * a*a*a*a*a*a
	   +0.00053485514521888073087240392846 * a*a*a*a*a*a*a
	   +0.00042891792986220150647633418796 * a*a*a*a*a*a*a*a
	   +0.00031639765900613633260381972850 * a*a*a*a*a*a*a*a*a
	   +0.00021908147678699592975840749194 * a*a*a*a*a*a*a*a*a*a
	   ;
  }
  else {
    double t   = 1./sqrt((z-1)*(z+1));
    double t2  = t*t;
    double rz  = sqrt(abs_zeta);
    double z32 = rz*rz*rz;
    double z92 = z32*z32*z32;
    double term1 = -t2*t * (30375. + 369603.*t2 + 765765.*t2*t2 + 425425.*t2*t2*t2)/414720.;
    double term2 = 85085./(663552.*z92);
    double term3 = -385./110592.*t*(3.+5.*t2)/(abs_zeta*abs_zeta*abs_zeta);
    double term4 = 5./55296.*t2*(81. + 462.*t2 + 385.*t2*t2)/z32;
    return (term1 + term2 + term3 + term4)/rz;
  }
}

/* checked OK [GJ] Thu Apr 30 22:27:36 MDT 1998 */
static double olver_B2(double z, double abs_zeta)
{
  if(z < 0.8) {
    double x = 5.*z/2. - 1.;
    return gsl_sf_cheb_eval(x,&B2_lt1_cs) / z;
  }
  else if(z <= 1.2) {
    double a = 1.-z;
    return   0.00055221307672129279005986982501
           + 0.00089586516310476929281129228969 * a
	   + 0.00067015003441569770883539158863 * a*a
	   + 0.00010166263361949045682945811828 * a*a*a
	   - 0.00044086345133806887291336488582 * a*a*a*a
	   - 0.00073963081508788743392883072523 * a*a*a*a*a
	   - 0.00076745494377839561259903887331 * a*a*a*a*a*a
	   - 0.00060829038106040362291568012663 * a*a*a*a*a*a*a
	   - 0.00037128707528893496121336168683 * a*a*a*a*a*a*a*a
	   - 0.00014116325105702609866850307176 * a*a*a*a*a*a*a*a*a
	   ;
  }
  else {
    double zi = 1./z;
    double x  = 12./5. * zi - 1.;
    return gsl_sf_cheb_eval(x,&B2_gt1_cs) * zi*zi*zi;
  }
}

/* checked OK [GJ] Thu Apr 30 21:01:51 MDT 1998  */
static double olver_B3(double z, double abs_zeta)
{
  if(z < 0.8) {
    double x = 5.*z/2. - 1.;
    return gsl_sf_cheb_eval(x,&B3_lt1_cs);
  }
  else if(z < 1.2) {
    double a = 1.-z;
    return -0.00047461779655995980754441833105
           -0.00095572913429464297452176811898 * a
	   -0.00080369634512082892655558133973 * a*a
	   -0.00000727921669154784138080600339 * a*a*a
	   +0.00093162500331581345235746518994 * a*a*a*a
	   +0.00149848796913751497227188612403 * a*a*a*a*a
	   +0.00148406039675949727870390426462 * a*a*a*a*a*a
	   ;
  }
  else {
    double x   = 12./(5.*z) - 1.;
    double zi2 = 1./(z*z);
    return gsl_sf_cheb_eval(x,&B3_gt1_cs) * zi2*zi2*zi2;
  }
}

/* checked OK [GJ] Thu Apr 30 20:08:13 MDT 1998 */
static double olver_A1(double z, double abs_zeta)
{
  if(z < 0.99) {
    double t = 1./sqrt((1-z)*(1+z));
    double rz = sqrt(abs_zeta);
    double t2 = t*t;
    double term1 =  t2*(81. - 462.*t2 + 385.*t2*t2)/1152.;
    double term2 = -455./(4608.*abs_zeta*abs_zeta*abs_zeta);
    double term3 =  7.*t*(-3 + 5.*t2)/(1152.*rz*rz*rz);
    return term1 + term2 + term3;
  }
  else if(z < 1.01) {
    double a = 1.-z;
    return -0.0044444444444444444444444444444
           -0.00184415584415584415584415584416  * a 
           +0.00056812076812076812076812076812  * a*a
	   +0.00168137865661675185484709294233  * a*a*a
	   +0.00186744042139000122193399504324  * a*a*a*a
	   +0.00161330105833747826430066790326  * a*a*a*a*a
	   +0.00123177312220625816558607537838  * a*a*a*a*a*a
	   +0.00087334711007377573881689318421  * a*a*a*a*a*a*a
	   +0.00059004942455353250141217015410  * a*a*a*a*a*a*a*a
	   ;
  }
  else {
    double t = 1./sqrt((z-1)*(z+1));
    double rz = sqrt(abs_zeta);
    double t2 = t*t;
    double term1 = -t2*(81. + 462.*t2 + 385.*t2*t2)/1152.;
    double term2 =  455./(4608.*abs_zeta*abs_zeta*abs_zeta);
    double term3 = -7.*t*( 3 + 5.*t2)/(1152.*rz*rz*rz);
    return term1 + term2 + term3;
  }
}

/* checked OK [GJ] Thu Apr 30 20:03:06 MDT 1998 */
static double olver_A2(double z, double abs_zeta)
{
  if(z < 0.88) {
    double t  = 1./sqrt((1-z)*(1+z));
    double t2 = t*t;
    double t4 = t2*t2;
    double t6 = t4*t2;
    double t8 = t4*t4;
    double rz = sqrt(abs_zeta);
    double z3 = abs_zeta*abs_zeta*abs_zeta;
    double z32 = rz*rz*rz;
    double z92 = z3*z32;
    double term1 = t4*(4465125. - 94121676.*t2 + 349922430.*t4 - 446185740.*t6  + 185910725.*t8)/39813120.;
    double term2 = -40415375./(127401984.*z3*z3);
    double term3 = -95095./15925248.*t*(3.-5.*t2)/z92;
    double term4 = -455./5308416. *t2*(81. - 462.*t2 + 385.*t4)/z3;
    double term5 = -7./19906560.*t*t2*(30375. - 369603.*t2  + 765765.*t4  - 425425.*t6)/z32;
    return term1 + term2 + term3 + term4 + term5;
  }
  else if(z < 1.12) {
    double a = 1.-z;
    return  0.000693735541354588973636592684210
           +0.000464483490365843307019777608010 * a
           -0.000289036254605598132482570468291 * a*a
	   -0.000874764943953712638574497548110 * a*a*a
	   -0.00102971637613986562996858467935  * a*a*a*a
	   -0.00083685732971381060058471403165  * a*a*a*a*a
	   -0.00048891089352721895499827012454  * a*a*a*a*a*a
	   -0.000144236747940817220502256810151 * a*a*a*a*a*a*a
	   +0.000114363800986163478038576460325 * a*a*a*a*a*a*a*a
	   +0.000266806881492777536223944807117 * a*a*a*a*a*a*a*a*a
	   -0.011975517576151069627471048587000 * a*a*a*a*a*a*a*a*a*a
	   ;
  }
  else {
    double t  = 1./sqrt((z-1)*(z+1));
    double t2 = t*t;
    double t4 = t2*t2;
    double t6 = t4*t2;
    double t8 = t4*t4;
    double rz = sqrt(abs_zeta);
    double z3 = abs_zeta*abs_zeta*abs_zeta;
    double z32 = rz*rz*rz;
    double z92 = z3*z32;
    double term1 = t4*(4465125. + 94121676.*t2 + 349922430.*t4 + 446185740.*t6  + 185910725.*t8)/39813120.;
    double term2 = -40415375./(127401984.*z3*z3);
    double term3 = +95095./15925248.*t*(3.+5.*t2)/z92;
    double term4 = -455./5308416. *t2*(81. + 462.*t2 + 385.*t4)/z3;
    double term5 = +7./19906560.*t*t2*(30375. + 369603.*t2  + 765765.*t4  + 425425.*t6)/z32;
    return term1 + term2 + term3 + term4 + term5;
  }
}

/* checked OK [GJ] Thu Apr 30 17:35:03 MDT 1998 */
static double olver_A3(double z, double abs_zeta)
{
  if(z < 0.9) {
    double x = 20.*z/9. - 1.;
    return gsl_sf_cheb_eval(x, &A3_lt1_cs);
  }
  else if(z < 1.1) {
    double a = 1.-z;
    return -0.00035421197145774384077112575920
           -0.000312322527890318832782774881353 * a
	   +0.000277947465383133980329617631915 * a*a
	   +0.000919803044747966977054155192400 * a*a*a
	   +0.001147600388275977640983696906320 * a*a*a*a
	   +0.00086923932612362574293177204454  * a*a*a*a*a
	   +0.000287392257282507334785281718027 * a*a*a*a*a*a
	   ;
  }
  else {
    double x   = 11./(5.*z) - 1.;
    double zi2 = 1./(z*z);
    return gsl_sf_cheb_eval(x, &A3_gt1_cs) * zi2*zi2*zi2;
  }
}

/* checked OK [GJ] Thu Apr 30 19:58:55 MDT 1998 */
static double olver_A4(double z, double abs_zeta)
{
  if(z < 0.8) {
    double x = 5.*z/2. - 1.;
    return gsl_sf_cheb_eval(x, &A4_lt1_cs);
  }
  else if(z < 1.2) {
    double a = 1.-z;
    return  0.00037819419920177291402661228437
           +0.00040494390552363233477213857527 * a
	   -0.00045764735528936113047289344569 * a*a
	   -0.00165361044229650225813161341879 * a*a*a
	   -0.00217527517983360049717137015539 * a*a*a*a
	   -0.00152003287866490735107772795537 * a*a*a*a*a
	   ;
  }
  else {
    double x   = 12./(5.*z) - 1.;
    double zi2 = 1./(z*z);
    return gsl_sf_cheb_eval(x, &A4_gt1_cs) * zi2*zi2*zi2*zi2;
  }
}


static inline double olver_Asum(double nu, double z, double abs_zeta)
{
  double nu2 = nu*nu;
  double A1 = olver_A1(z, abs_zeta);
  double A2 = olver_A2(z, abs_zeta);
  double A3 = olver_A3(z, abs_zeta);
  double A4 = olver_A4(z, abs_zeta);
  return 1. + A1/nu2 + A2/(nu2*nu2) + A3/(nu2*nu2*nu2) + A4/(nu2*nu2*nu2*nu2);
}

static inline double olver_Bsum(double nu, double z, double abs_zeta)
{
  double nu2 = nu*nu;
  double B0 = olver_B0(z, abs_zeta);
  double B1 = olver_B1(z, abs_zeta);
  double B2 = olver_B2(z, abs_zeta);
  double B3 = olver_B3(z, abs_zeta);
  return B0 + B1/nu2 + B2/(nu2*nu2) + B3/(nu2*nu2*nu2*nu2);
}


/* uniform asymptotic, nu -> Inf, [Abramowitz+Stegun, 9.3.35]
 *
 * error:
 *    nu =  2: uniformly good to >  6D
 *    nu =  5: uniformly good to >  8D
 *    nu = 10: uniformly good to > 10D
 *    nu = 20: uniformly good to > 13D
 *
 * checked OK [GJ] Sun May  3 22:36:29 EDT 1998 
 */
int gsl_sf_bessel_Jnu_asymp_Olver_impl(double nu, double x, double * result)
{
  if(x <= 0. || nu <= 0.) {
    return GSL_EDOM;
  }  
  else {
    double zeta, abs_zeta;
    double arg;
    double pre;
    double asum, bsum;
    double ai, aip;
    double z = x/nu;
    double crnu = pow(nu, 1./3.);
    double rt   = sqrt(fabs(1.-z)*(1+z));
    
    if(fabs(1-z) < GSL_ROOT4_MACH_EPS) {
      /* z near 1 */
      double a = 1.-z;
      pre  =  1.25992104989487316476721060728
             +0.37797631496846194943016318218  * a
	     +0.230385563409348235843147082474 * a*a
	     +0.165909603649648694839821892031 * a*a*a
	     ;
      zeta = a * pre;
      pre  = sqrt(sqrt(4.*pre/(1+z)));
      abs_zeta = fabs(zeta);
    }
    else if(z < 1.) {
      /* z < 1 */
      abs_zeta = pow(1.5*(log((1+rt)/z) - rt), 2./3.);
      zeta = abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }
    else {
      /* z > 1 */
      abs_zeta = pow(1.5*(rt - acos(1./z)), 2./3.);
      zeta = -abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }

    asum = olver_Asum(nu, z, abs_zeta);
    bsum = olver_Bsum(nu, z, abs_zeta);
    arg  = crnu*crnu * zeta;
    gsl_sf_airy_Ai_impl(arg, &ai);
    gsl_sf_airy_Ai_deriv_impl(arg, &aip);

    *result = pre * (ai*asum/crnu + aip*bsum/(nu*crnu*crnu));
    return GSL_SUCCESS;
  }
}

/* uniform asymptotic, nu -> Inf,  [Abramowitz+Stegun, 9.3.36]
 *
 * error:
 *    nu =  2: uniformly good to >  6D
 *    nu =  5: uniformly good to >  8D
 *    nu = 10: uniformly good to > 10D
 *    nu = 20: uniformly good to > 13D
 *
 * checked OK [GJ] Sun May  3 22:59:22 EDT 1998 
 */
int gsl_sf_bessel_Ynu_asymp_Olver_impl(double nu, double x, double * result)
{
  if(x <= 0. || nu <= 0.) {
    return GSL_EDOM;
  }  
  else {
    double zeta, abs_zeta;
    double arg;
    double pre;
    double asum, bsum;
    double bi, bip;
    double z = x/nu;
    double crnu = pow(nu, 1./3.);
    double rt   = sqrt(fabs(1.-z)*(1+z));
    
    if(fabs(1-z) < GSL_ROOT4_MACH_EPS) {
      /* z near 1 */
      double a = 1.-z;
      pre  =  1.25992104989487316476721060728
             +0.37797631496846194943016318218  * a
	     +0.230385563409348235843147082474 * a*a
	     +0.165909603649648694839821892031 * a*a*a
	     ;
      zeta = a * pre;
      pre  = sqrt(sqrt(4.*pre/(1+z)));
      abs_zeta = fabs(zeta);
    }
    else if(z < 1.) {
      /* z < 1 */
      abs_zeta = pow(1.5*(log((1+rt)/z) - rt), 2./3.);
      zeta = abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }
    else {
      /* z > 1 */
      abs_zeta = pow(1.5*(rt - acos(1./z)), 2./3.);
      zeta = -abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }

    asum = olver_Asum(nu, z, abs_zeta);
    bsum = olver_Bsum(nu, z, abs_zeta);
    arg  = crnu*crnu * zeta;
    gsl_sf_airy_Bi_impl(arg, &bi);
    gsl_sf_airy_Bi_deriv_impl(arg, &bip);

    *result = -pre * (bi*asum/crnu + bip*bsum/(nu*crnu*crnu));
    return GSL_SUCCESS;
  }
}
