! molecules_module.f90
module molecules_module
  use amr_parameters, only: dp
  use safe_math, only: safe_exp
  use rtz_module, only: G0_draine
  implicit none

  private  ! everything is private by default
  public :: alpha_H2, beta_H2_umist, beta_H2, beta_H2_krome, alpha_CO, beta_CO, alpha_H2_prim, alpha_H2_dust
  public :: comp_SH2, comp_Sd, comp_SCO, initialize_SCO_table, lw_transmission, tau_dust_lw_mw, co_dust_to_lw
  public :: fuv_dust_to_lw

  real(dp), dimension(52, 2):: sco_table ! self shielding coefficient vs. CO column density
  real(dp), dimension(43, 2):: sh2_table ! self shielding coefficient vs. H2 column density

  ! Dust attenuation in the Lyman-Werner band (1000 A) per H nucleus, for MW dust (R_V = 3.1), on the
  ! total hydrogen column N(H) + 2N(H2) + N(H+): Draine & Bertoldi (1996, Sect. 2.4), an effective
  ! attenuation cross-section (~0.76 of the extinction one, as forward-scattered light is not removed)
  real(dp), parameter :: sigma_dust_lw_mw = 2.0d-21
  ! CO and H2 dust attenuation in units of A_V: CO exp(-3.53 A_V) (Visser et al. 2009; van Dishoeck
  ! et al. 2006), H2 sigma_dust_lw_mw N_H = 3.74 A_V for N_H/A_V = 1.87e21 cm^-2 mag^-1
  real(dp), parameter :: co_dust_to_lw = 3.53d0 / 3.74d0
  ! The FUV (6-13.6 eV) field's dust attenuation, exp(-1.8 A_V) (Tielens & Hollenbach 1985), for MW dust
  real(dp), parameter :: fuv_dust_to_lw = 1.8d0 / 3.74d0

  ! e^x E2(x) for 1 < x < 1024 (e2_scaled): Chebyshev coefficients on the sub-intervals k = 4 e + s,
  ! [2^e (1 + s/4), 2^e (1 + (s+1)/4)], e = 0..9, s = 0..3 (fitted to mpmath; relative error < 4e-16)
  integer, parameter :: E2_NC = 15, E2_NI = 40
  real(dp), parameter :: E2_CHEB(0:E2_NC-1, 0:E2_NI-1) = reshape([ &
       3.81763510007170681d-01, -2.11987461473099013d-02, 6.66392660023101303d-04, &
       -2.30907959755498501d-05, 8.62459451977421257d-07, -3.41306824135789781d-08, &
       1.41279795793731402d-09, -6.05938506084286222d-11, 2.67387520633926742d-12, &
       -1.20762527308868542d-13, 5.55985005241087548d-15, -2.60129900214044407d-16, &
       1.23384164350579036d-17, -5.92148869185377367d-19, 2.87095879302252799d-20, &
       3.43957231796447049d-01, -1.67835831000091515d-02, 4.54294819733923262d-04, &
       -1.33527481887875754d-05, 4.18644246457620606d-07, -1.38073846097764668d-08, &
       4.73993422803296969d-10, -1.68021347708667901d-11, 6.11323298981706285d-13, &
       -2.27245828831130140d-14, 8.59999681368906907d-16, -3.30424293637186882d-17, &
       1.28605696926891891d-18, -5.06167327646441715d-20, 2.01163027454355636d-21, &
       3.13585689353303720d-01, -1.36939171454812780d-02, 3.26782359098568877d-04, &
       -8.37352105068389131d-06, 2.27048646469836361d-07, -6.43976789461903433d-09, &
       1.89360008028624518d-10, -5.73332065033383445d-12, 1.77806064272521659d-13, &
       -5.62531896104191078d-15, 1.80979794701758781d-16, -5.90615953873498477d-18, &
       1.95118924219256405d-19, -6.51486546230452338d-21, 2.19555294947528783d-22, &
       2.88531929872702819d-01, -1.14282452964928876d-02, 2.44520698610506165d-04, &
       -5.56918286643584194d-06, 1.33373978082364994d-07, -3.32586724843230762d-09, &
       8.56985395418428527d-11, -2.26829581544321721d-12, 6.13873960763445001d-14, &
       -1.69255252435717934d-15, 4.74075601615350150d-17, -1.34587372780576648d-18, &
       3.86555528924284889d-20, -1.12154612169485935d-21, 3.28308070495305272d-23, &
       2.58590925533045912d-01, -1.80522575926450601d-02, 6.72079491946174024d-04, &
       -2.63751951076197584d-05, 1.08047336361365198d-06, -4.58417811775942225d-08, &
       2.00184719211521680d-09, -8.95336503580305881d-11, 4.08544178244309012d-12, &
       -1.89603654034017682d-13, 8.92757702816521401d-15, -4.25627766598427618d-16, &
       2.05128422171450655d-17, -9.98013521573811332d-19, 4.89638456486404444d-20, &
       2.27030920401445013d-01, -1.37038365316549948d-02, 4.35675602757850907d-04, &
       -1.44606005361979867d-05, 4.97364671329190625d-07, -1.76195549552607011d-08, &
       6.39771781037569372d-10, -2.37171478813067438d-11, 8.94823916722824230d-13, &
       -3.42723389883742323d-14, 1.32978262232379648d-15, -5.21807535981441443d-17, &
       2.06786924129251556d-18, -8.26629598474974782d-20, 3.33002592835383767d-21, &
       2.02639835255185169d-01, -1.07995451374009264d-02, 3.00466154111953823d-04, &
       -8.66691101351896439d-06, 2.57666157583366930d-07, -7.85735042597861928d-09, &
       2.44798925323852988d-10, -7.76721102116653259d-12, 2.50326975115522281d-13, &
       -8.17725541124116053d-15, 2.70270879196427767d-16, -9.02501714118190495d-18, &
       3.04105340824743020d-19, -1.03295598559854561d-20, 3.53381818106796646d-22, &
       1.83157846245583233d-01, -8.75192285629800872d-03, 2.16878149375730744d-04, &
       -5.54274673985613996d-06, 1.45399746531940764d-07, -3.89959274709054545d-09, &
       1.06581354715287815d-10, -2.96066135573518084d-12, 8.34035465070757467d-14, &
       -2.37835204132127625d-15, 6.85494809741836906d-17, -1.99442014720541292d-18, &
       5.85125061825687211d-20, -1.72944316479810983d-21, 5.14580577511127314d-23, &
       1.60681937349291010d-01, -1.33390953789767168d-02, 5.70200320725052401d-04, &
       -2.49943456005639935d-05, 1.11944301271170651d-06, -5.10735872051535655d-08, &
       2.36772105748868166d-09, -1.11299082951921936d-10, 5.29561388497682581d-12, &
       -2.54662151022998433d-13, 1.23622128956981849d-14, -6.05138320478750330d-16, &
       2.98435477077286135d-17, -1.48166593652493186d-18, 7.40060560143489304d-20, &
       1.37789407333192226d-01, -9.73446305433203242d-03, 3.51777425315460263d-04, &
       -1.29654078435975294d-05, 4.86113615416961220d-07, -1.84997915127618632d-08, &
       7.13287911525970718d-10, -2.78193889957896411d-11, 1.09606045019235073d-12, &
       -4.35744146582180590d-14, 1.74629094455370949d-15, -7.04899755878585718d-17, &
       2.86385841713322420d-18, -1.17035790757123955d-19, 4.80832153461712195d-21, &
       1.20719958116962733d-01, -7.43351496242761374d-03, 2.33077734089618592d-04, &
       -7.42574759735946181d-06, 2.39919255941795861d-07, -7.84785317534318992d-09, &
       2.59524128362706347d-10, -8.66592744773847099d-12, 2.91882004016759426d-13, &
       -9.90746238388514295d-15, 3.38642621958526799d-16, -1.16479925924292389d-17, &
       4.02935042207308721d-19, -1.40110235319463921d-20, 4.89507434430227558d-22, &
       1.07474746402021148d-01, -5.87014575853084079d-03, 1.62730282267794135d-04, &
       -4.57125263894638240d-06, 1.29925238014479218d-07, -3.73144822027630387d-09, &
       1.08167918618180648d-10, -3.16179982851338141d-12, 9.31148887110654219d-14, &
       -2.76079818733211753d-15, 8.23575908000754738d-17, -2.47049176940726349d-18, &
       7.44836528585241711d-20, -2.25604938148519361d-21, 6.86246810306511731d-23, &
       9.26526352896867439d-02, -8.67838087439816268d-03, 4.11184480745827681d-04, &
       -1.96853464945509828d-05, 9.51235365274286693d-07, -4.63512810946163094d-08, &
       2.27565545509153632d-09, -1.12488791537585156d-10, 5.59494062968077350d-12, &
       -2.79849081064891367d-13, 1.40695984282984302d-14, -7.10690969202048794d-16, &
       3.60540700421478543d-17, -1.83634273453244671d-18, 9.38742831481667982d-20, &
       7.79830912069654214d-02, -6.13040180405167046d-03, 2.43048820912088219d-04, &
       -9.71147697622520533d-06, 3.90793817086345318d-07, -1.58270949423518481d-08, &
       6.44760877689411730d-10, -2.64071230926415071d-11, 1.08686093512281743d-12, &
       -4.49349502370836721d-14, 1.86550972240184269d-15, -7.77455595155196249d-17, &
       3.25158121856384478d-18, -1.36440455864895678d-19, 5.74276233833974260d-21, &
       6.73596345026768878d-02, -4.56547587750313726d-03, 1.55754627059153029d-04, &
       -5.34616332724345697d-06, 1.84529246566764055d-07, -6.40190988868415127d-09, &
       2.23149680043361657d-10, -7.81205203133559399d-12, 2.74582236595965864d-13, &
       -9.68699443148711724d-15, 3.42924589954432013d-16, -1.21785658450635356d-17, &
       4.33797043196858922d-19, -1.54947215553860217d-20, 5.54893201676196843d-22, &
       5.93015822477188786d-02, -3.53401176783157750d-03, 1.05865140096617728d-04, &
       -3.18688740703388522d-06, 9.63703635880548887d-08, -2.92640370394817327d-09, &
       8.92081693580590652d-11, -2.72918787059550678d-12, 8.37743412113429176d-14, &
       -2.57951502367315988d-15, 7.96570790043143272d-17, -2.46653873794896743d-18, &
       7.65691793226039629d-20, -2.38261713552398924d-21, 7.43064101831488640d-23, &
       5.04896823542088785d-02, -5.10737210256083287d-03, 2.59355924219295943d-04, &
       -1.32192051336137409d-05, 6.76103367621253043d-07, -3.46910729482209979d-08, &
       1.78536350317493902d-09, -9.21415344296104664d-11, 4.76789843558768022d-12, &
       -2.47326576738156791d-13, 1.28594436558619945d-14, -6.70072584577760710d-16, &
       3.49876882856982507d-17, -1.83042585887202401d-18, 9.59370129920480795d-20, &
       4.19499622612159243d-02, -3.52412994308717734d-03, 1.48450516603297199d-04, &
       -6.27009431921501650d-06, 2.65496268463785376d-07, -1.12685982673007257d-08, &
       4.79345668162279664d-10, -2.04333767208839404d-11, 8.72756577912207033d-13, &
       -3.73474411223252211d-14, 1.60103150260059045d-15, -6.87495662605812744d-17, &
       2.95688038383410294d-18, -1.27366891069858932d-19, 5.49421992699691641d-21, &
       3.58928071822226455d-02, -2.57912233810678529d-03, 9.28619514818470915d-05, &
       -3.35030028115616491d-06, 1.21105047503008302d-07, -4.38560302052739882d-09, &
       1.59090443151850528d-10, -5.78054242287811745d-12, 2.10362094698707370d-13, &
       -7.66668917547211362d-15, 2.79807274028164092d-16, -1.02256812451764543d-17, &
       3.74180203377291923d-19, -1.37088045723736218d-20, 5.02834520680428175d-22, &
       3.13697144003640890d-02, -1.96964591928666926d-03, 6.19386466104563971d-05, &
       -1.95085625805971812d-06, 6.15382747534657104d-08, -1.94397443015209281d-09, &
       6.14938538795166339d-11, -1.94778687224793853d-12, 6.17723437269847756d-14, &
       -1.96139746200819369d-15, 6.23495700930739964d-17, -1.98415889270321243d-18, &
       6.32083758412237221d-20, -2.01562275230721514d-21, 6.43372707290565646d-23, &
       2.64987528583138891d-02, -2.80445588392760469d-03, 1.48583804555639488d-04, &
       -7.88135980029304578d-06, 4.18519935242344908d-07, -2.22482684695984255d-08, &
       1.18391979483681877d-09, -6.30632514177598613d-11, 3.36233877263528499d-12, &
       -1.79432693158032379d-13, 9.58388331505257052d-15, -5.12326651335155460d-16, &
       2.74097001440771303d-17, -1.46757617217363494d-18, 7.86364395600301831d-20, &
       2.18419919839617782d-02, -1.90633710040493841d-03, 8.32616941569437597d-05, &
       -3.63953035875760983d-06, 1.59216172332372158d-07, -6.97041378968564102d-09, &
       3.05385151712827976d-10, -1.33888799475604212d-11, 5.87404249261845641d-13, &
       -2.57879287821756756d-14, 1.13285125627146119d-15, -4.97963283509540748d-17, &
       2.19018664945913454d-18, -9.63864692491292243d-20, 4.24418486085708489d-21, &
       1.85820124797359716d-02, -1.38015024589927633d-03, 5.12862297449022508d-05, &
       -1.90694118702867939d-06, 7.09459601046093293d-08, -2.64097313173326755d-09, &
       9.83645658423551265d-11, -3.66559543805364661d-12, 1.36670448819297337d-13, &
       -5.09826402882380780d-15, 1.90274807874116649d-16, -7.10469510527696812d-18, &
       2.65404907256276789d-19, -9.91896152297995323d-21, 3.70861579467420333d-22, &
       1.61709505129545260d-02, -1.04542059831232723d-03, 3.38083657072489743d-05, &
       -1.09385317959093350d-06, 3.54070686815488772d-08, -1.14659904069983162d-09, &
       3.71465620668227598d-11, -1.20394354159789944d-12, 3.90363538743562464d-14, &
       -1.26620309158994854d-15, 4.10869623719760634d-17, -1.33372702108137168d-18, &
       4.33099800993539601d-20, -1.40689975917731352d-21, 4.57181476186410358d-23, &
       1.35981632173124214d-02, -1.47544602158956721d-03, 8.00728284114231472d-05, &
       -4.34702524733772289d-06, 2.36070009241505374d-07, -1.28241256459810759d-08, &
       6.96866762535605571d-10, -3.78794472588889977d-11, 2.05961654350483010d-12, &
       -1.12019860652352677d-13, 6.09434113647636042d-15, -3.31649097112301775d-16, &
       1.80529651092465283d-17, -9.82954362551543822d-19, 5.35341312282996949d-20, &
       1.11580188404434271d-02, -9.94271569034954053d-04, 4.43092542410824161d-05, &
       -1.97507359486373186d-06, 8.80581478243187238d-08, -3.92691307526119157d-09, &
       1.75156645811176078d-10, -7.81436529134464992d-12, 3.48699057175586660d-13, &
       -1.55630990262966714d-14, 6.94749029215533914d-16, -3.10202104190712646d-17, &
       1.38530327564028650d-18, -6.18766459846044804d-20, 2.76432399269982816d-21, &
       9.46262283204899433d-03, -7.15429860957970021d-04, 2.70499375268229001d-05, &
       -1.02291100144885819d-06, 3.86883811914453754d-08, -1.46350149442677635d-09, &
       5.53699961722714611d-11, -2.09518977274840746d-12, 7.92936943360604252d-14, &
       -3.00136799892617854d-15, 1.13622429730125179d-16, -4.30201650184623606d-18, &
       1.62907927545756300d-19, -6.16983441331297197d-21, 2.33703362363381429d-22, &
       8.21549913530385256d-03, -5.39444017443323947d-04, 1.77126982509561215d-05, &
       -5.81671952235217218d-07, 1.91040629355361411d-08, -6.27518822478413582d-10, &
       2.06148591978002658d-11, -6.77307374733810682d-13, 2.22557526185972246d-14, &
       -7.31390080650883043d-16, 2.40383919825838101d-17, -7.90151934035031518d-19, &
       2.59755018713965291d-20, -8.54013516101016392d-22, 2.80809841400878987d-23, &
       6.89147595150569128d-03, -7.57653160771419565d-04, 4.16522309851241216d-05, &
       -2.29005130608511037d-06, 1.25918855345618567d-07, -6.92427706357419380d-09, &
       3.80798944858993802d-10, -2.09437375871239321d-11, 1.15199162911531862d-12, &
       -6.33695588050158960d-14, 3.48616402677933168d-15, -1.91800741378476397d-16, &
       1.05532843936249879d-17, -5.80710222277647964d-19, 3.19569519608315450d-20, &
       5.64117137129405404d-03, -5.08163807139600334d-04, 2.28894267136897368d-05, &
       -1.03108046262143020d-06, 4.64490005895889309d-08, -2.09259935000504559d-09, &
       9.42804007651677428d-11, -4.24797650345601456d-12, 1.91411434667270400d-13, &
       -8.62538473752906172d-15, 3.88699195038684292d-16, -1.75175401122999812d-17, &
       7.89508235317721732d-19, -3.55847588342810294d-20, 1.60396528903276329d-21, &
       4.77600433497256047d-03, -3.64448527997177953d-04, 1.39058313550154349d-05, &
       -5.30611678601725598d-07, 2.02476938391379171d-08, -7.72668253822995052d-10, &
       2.94869019712778995d-11, -1.12533967090398149d-12, 4.29493260334681031d-14, &
       -1.63925741052663586d-15, 6.25685190477475191d-17, -2.38826412337155133d-18, &
       9.11646431314317425d-20, -3.48006988918057952d-21, 1.32851600549993622d-22, &
       4.14143452737674757d-03, -2.74132346294764142d-04, 9.07307100889790180d-06, &
       -3.00305135616585936d-07, 9.93998014291638866d-09, -3.29020129540368952d-10, &
       1.08911440903031817d-11, -3.60527594119421824d-13, 1.19348613855832843d-14, &
       -3.95102704456659782d-16, 1.30802568653524666d-17, -4.33048046833881438d-19, &
       1.43373647008530443d-20, -4.74696374168547854d-22, 1.57172198398213874d-23, &
       3.46954708809307354d-03, -3.84037492979131146d-04, 2.12546869826425054d-05, &
       -1.17637551213329014d-06, 6.51099472736635719d-08, -3.60378398033781982d-09, &
       1.99471167163920410d-10, -1.10410763422678524d-11, 6.11156642766561793d-13, &
       -3.38301119723947070d-14, 1.87268203629616150d-15, -1.03665512898996230d-16, &
       5.73870817880172758d-18, -3.17689969025143733d-19, 1.75874276280927719d-20, &
       2.83651736397884640d-03, -2.56942525228411945d-04, 1.16375994663952422d-05, &
       -5.27105633392465093d-07, 2.38747430249700759d-08, -1.08140036705622902d-09, &
       4.89825106167410921d-11, -2.21871875686366295d-12, 1.00500938412710009d-13, &
       -4.55244378053076730d-15, 2.06217559081492242d-16, -9.34142719666373333d-18, &
       4.23162652089936837d-19, -1.91693725400065556d-20, 8.68390173625690541d-22, &
       2.39941357392347079d-03, -1.83961303057813229d-04, 7.05217006908743896d-06, &
       -2.70348512169584186d-07, 1.03640637156002142d-08, -3.97320547787176623d-10, &
       1.52319959330616665d-11, -5.83952363155740869d-13, 2.23873566922749841d-14, &
       -8.58287897062243407d-16, 3.29054555142049171d-17, -1.26155896167127246d-18, &
       4.83673066529263923d-20, -1.85438941342877699d-21, 7.10975502101843484d-23, &
       2.07929422879007172d-03, -1.38199257868465739d-04, 4.59271173139400624d-06, &
       -1.52628753491702560d-07, 5.07232656789864466d-09, -1.68570556255829183d-10, &
       5.60221642328713300d-12, -1.86183709625117290d-13, 6.18766931280244883d-15, &
       -2.05644035271423927d-16, 6.83453118431086821d-18, -2.27145908141090215d-19, &
       7.54926495276400804d-21, -2.50904219810363718d-22, 8.33901579741495363d-24, &
       1.74081892037828258d-03, -1.93351841276090701d-04, 1.07378080314124906d-05, &
       -5.96328484843058775d-07, 3.31175393740891112d-08, -1.83921773564191441d-09, &
       1.02143515410337421d-10, -5.67271566490477399d-12, 3.15045864513915837d-13, &
       -1.74968187885027779d-14, 9.71732972963508608d-16, -5.39681227123743511d-17, &
       2.99729991163498557d-18, -1.66466021325011195d-19, 9.24535303114911273d-21, &
       1.42229242556055828d-03, -1.29200046008730161d-04, 5.86824430517642689d-06, &
       -2.66535733122852053d-07, 1.21061047784970608d-08, -5.49863818225317351d-10, &
       2.49751198669913490d-11, -1.13438829601331116d-12, 5.15249532010909022d-14, &
       -2.34031967276926327d-15, 1.06300293662627066d-16, -4.82831331446841875d-18, &
       2.19309840102932446d-19, -9.96144822034190844d-21, 4.52468687959630484d-22, &
       1.20259040194530746d-03, -9.24219126698613862d-05, 3.55143133423712936d-06, &
       -1.36468725777274925d-07, 5.24401701482862874d-09, -2.01509852939420287d-10, &
       7.74336471897698791d-12, -2.97553030044197972d-13, 1.14340547761424769d-14, &
       -4.39377067885866468d-16, 1.68840149480437635d-17, -6.48806646874588768d-19, &
       2.49319398777211735d-20, -9.58071820374922378d-22, 3.68163961431845537d-23, &
       1.04181138306865635d-03, -6.93868588884038963d-05, 2.31066133990211741d-06, &
       -7.69478170287144630d-08, 2.56246086940850569d-09, -8.53334121541814477d-11, &
       2.84172410064324852d-12, -9.46336599278550923d-14, 3.15144903453421879d-15, &
       -1.04948409934665590d-16, 3.49496120656452771d-18, -1.16388420843565192d-19, &
       3.87594949809456106d-21, -1.29076553207179656d-22, 4.29850586607345343d-24 &
       ], [E2_NC, E2_NI])

CONTAINS

FUNCTION alpha_H2_prim(T, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) result(rate)
  implicit none

  real(dp), intent(in) :: T, xe, H2_cosmic_ray_ionization_rate, G0
  real(dp), intent(in) :: xHI, xHII, nH
  real(dp) :: rate
  real(dp) :: logT, lnTe, logT2, Te
  real(dp) :: k1, k2, k5, k13, k14, k15, k_hm_cr, k_hm_gamma
  real(dp) :: cr_photo_destruction, collisional_destruction

  ! H- channel for H2 formation
  rate = 0.d0

  ! Primordial channel
  logT = log10(T)
  Te = T*8.621738d-5 ! K -> eV
  lnTe = log(Te)

  ! Creation and destruction channels of H- included with updated rates from Glover et al. 2010
  ! H + e- -> H- + gamma
  if (T .lt. 6000.d0) then
     k1 = 10.d0**(-17.845d0 + logT * (0.762d0 + logT * (0.1523d0 - 0.03274d0 * logT)))
  else
     logT2 = logT * logT
     k1 = 10.d0**(-16.42d0 + logT2 * (0.1998d0 + logT2 * (-5.447d-3 + 4.0415d-5 * logT2)))
  end if

  ! H- + H -> H2 + e
  k2 = 4.0d-9*(max(T,300.d0)**(-0.17d0))

  ! H- + H+ -> H + H
  k5 = 2.4d-6/sqrt(T)*(1.d0 + T/20000.d0)

  ! H- + CR --> H + e-
  k_hm_cr = 1.28d-13 * (H2_cosmic_ray_ionization_rate / 1.d-16)

  ! H- + gamma --> H + e-
  k_hm_gamma = 5.9d-9 * G0

  ! H- + e -> H + e + e
  k13 = -1.801849334d1 + lnTe * (2.36085220d0 + lnTe * (-2.82744300d-1 + lnTe * (1.62331664d-2 + lnTe * (-3.36501203d-2 + lnTe * (1.17832978d-2 + lnTe * (-1.65619470d-3 + lnTe * (1.06827520d-4 - 2.63128581d-6 * lnTe)))))))
  k13 = safe_exp(max(k13,-92.d0)) ! Max needed to prevent result from diverging at low temperature

  ! H- + H --> H + H + e-
  ! I think this reaction was broken in glover so I took the results from
  ! https://www.aanda.org/articles/aa/pdf/2016/02/aa27262-15.pdf Table A1
  ! Harley added the fudge factor for continuity
  k14 = 1.357772745525155d0 * 2.5634d-15 * (Te**1.78186d0) ! Note that T must be in eV for this reaction to make sense
  if (T .gt. 1160.d0) then
     k14 = -3.388464953d1 + lnTe * (1.13944933d0 + lnTe * (-1.4210135d-1 + lnTe * (8.4644554d-3 + lnTe * (-1.4328641d-3 + lnTe * (2.0122503d-4 + lnTe * (8.6639632d-5 + lnTe * (-2.5850097d-5 + lnTe * (2.4555012d-6 - 8.0683825d-8 * lnTe))))))))
     k14 = safe_exp(k14)
  end if

  ! H- + H+ --> H2+ + e- --> H + H (via recombinative dissociation)
  k15 = 6.9d-9 * (T**(-0.35d0))
  if (T .gt. 8000.d0) then
     k15 = 9.6d-7 * (T**(-0.9d0))
  end if

  ! Correct steady state denominator terms scaled to cm^3 s^-1
  cr_photo_destruction = (k_hm_cr + k_hm_gamma) / (nH + 1d-40)
  collisional_destruction = k2 * xHI + k5 * xHII + k13 * xe + k14 * xHI + k15 * xHII

  rate = rate + k1 * k2 * xe * xHI / (collisional_destruction + cr_photo_destruction)

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_H2_prim

FUNCTION alpha_H2_dust(T, dust_to_gas_mass_ratio_over_mw) result(rate)
  ! Formation on dust
  use rt_parameters, only: rtz_H2_clumping
  implicit none

  real(dp), intent(in) :: T, dust_to_gas_mass_ratio_over_mw
  real(dp) :: rate
  real(dp) :: T2

  T2 = T / 100.d0
  rate = dust_to_gas_mass_ratio_over_mw * (3.5d-17) * rtz_H2_clumping * sqrt(min(T2,1.d2))

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_H2_dust

FUNCTION alpha_H2(T, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) result(rate)
  ! Creation rate of molecular hydrogen
  ! We consider both the primordial channel (via H-) as well
  ! as formation on dust
  implicit none
  real(dp), intent(in) :: T, dust_to_gas_mass_ratio_over_mw, xe
  real(dp), intent(in) :: H2_cosmic_ray_ionization_rate, G0, xHI
  real(dp), intent(in) :: xHII, nH
  real(dp) :: rate

  rate = 0.d0

  ! Formation rate on dust. Consider only HI
  if (dust_to_gas_mass_ratio_over_mw.gt.0d0) then
     rate = rate + alpha_H2_dust(T, dust_to_gas_mass_ratio_over_mw) * xHI * nH
  end if

  ! Primordial H- channel
  rate = rate + alpha_H2_prim(T, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) * xHI * nH

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_H2

FUNCTION beta_H2_umist(T, nH, ne, nH2) result(rate)
  ! H2 destruction from umist
  implicit none

  real(dp), intent(in) :: T, nH, ne, nH2
  real(dp) :: rate
  real(dp) :: T_loc

  rate = 0.d0

  !H2 + H2 --> H2 + H + H 
  T_loc = max(min(T,41000.d0),2803.d0)
  rate = rate + (1.00d-8 * ((T_loc/300d0)**0.0d0) * safe_exp(-84100.d0/T_loc) * nH2)

  !H2 + e- --> H + H + e-
  T_loc = max(min(T,41000.d0),3400.d0)
  rate = rate + (3.22d-9 * ((T_loc/300d0)**0.35d0) * safe_exp(-102000.d0/T_loc) * ne)

  !H2 + H --> H + H + H
  T_loc = max(min(T,41000.d0),1833.d0)
  rate = rate + (4.67d-7 * ((T_loc/300d0)**(-1.d0)) * safe_exp(-55000.d0/T_loc) * nH)

  rate = MAX(rate,1.d-100)

END FUNCTION beta_H2_umist

FUNCTION beta_H2_krome(T, nH, ne, nH2, nHe) result(rate)
  ! From bovino 2016
  ! https://www.aanda.org/articles/aa/pdf/2016/06/aa28158-16.pdf
  implicit none
  real(dp), intent(in) :: T, nH, ne, nH2, nHe
  real(dp) :: rate

  rate = 0.d0

  ! H2 + H --> H + H + H (k18)
  rate = rate + ((6.67d-12 * sqrt(T) * safe_exp(-1.d0 * (1.d0 + (63593.d0/T)))) * nH)

  ! H2 + H2 --> H2 + H + H (k19)
  rate = rate + ((5.996d-30 * (T**4.1881d0) * ((1.d0 + (6.761d-6 * T))**(-5.6881d0)) * safe_exp(-54657.4d0/T)) * nH2)

  ! H2 + e- --> H + H + e- (k22)
  rate = rate + (4.38d-10 * (T**0.35d0) * safe_exp(-102000.d0/T) * ne)

  ! H2 + He --> H + H + He
  rate = rate + ((10.d0**(-27.029d0 + (3.801d0 * log10(T)) - (29487.d0/T))) * nHe)

  rate = MAX(rate,1.d-100)

END FUNCTION beta_H2_krome

FUNCTION beta_H2(T, nH, xHI, xH2, xHe, ne, nHI, nH2, nHeI) result(rate)
  ! Returns the collisional dissociation rates of H2 for four different
  ! reactions [cm3s-1] from Glover & Abel (2008)
  ! http://mnras.oxfordjournals.org/content/388/4/1627.full.pdf 
  implicit none

  real(dp), intent(in):: T, nH, xHI, xH2, xHe
  real(dp), intent(in):: ne, nHI, nH2, nHeI
  real(dp):: rate
  real(dp):: T4, ncrH, ncrH2, ncrHe, invncr
  real(dp):: LTEfac, NLTEfac
  real(dp):: k8, k9, k9L, k10, k10L, k11, k11L
  real(dp):: lk8, lk9, lk10, lk11

  T4 = T / 1d4

  ! Critical number densities.  See eqns 15, 16, and 17
  ncrH =  10.d0**(3.d0 - 0.416d0*log10(T4) - 0.327d0*log10(T4)*log10(T4))
  ncrH2 = 10.d0**(4.845d0 - 1.3d0*log10(T4) + 1.62d0*log10(T4)*log10(T4))
  ncrHe = 10.d0**(5.0792d0*(1.d0 - 1.23d-5*(T - 2000.d0)))

  ! 1/ncr.  see eqn 14
  invncr = (xHI/ncrH) + (xH2/ncrH2) + (xHe/ncrHe)

  ! prefactors for LTE and NLTE collision rates  see eqn 13
  LTEfac = (nH*invncr)/(1.d0 + (nH*invncr))
  NLTEfac = 1.d0/(1.d0 + (nH*invncr))
  if (ne .gt. nHI) then
     LTEfac = 1.d0
     NLTEfac = 0.d0
  end if

  !reaction rates from the appendix:
  !H2 + e- --> H + H + e-
  k8 = 3.73d-9 * (T**0.1121d0) * safe_exp(-99430.d0/T) ! Glover et al. (2010)

  !H2 + H --> H + H + H
  k9 = (6.67d-12)*sqrt(T)*safe_exp(-1.d0*(1.d0 + (63593.d0/T)))
  k9L = (3.52d-9)*safe_exp(-43900.d0/T)

  !H2 + H2 --> H2 + H + H
  k10 = ( (5.996d-30*(T**4.1881d0)) / ((1.d0 + 6.761d-6*T)**5.6881d0)) * safe_exp(-54657.4d0/T)
  k10L = (1.3d-9)*safe_exp(-53300.d0/T)

  !H2 + He --> H + H + He
  k11 = 10.d0**(-27.029d0 + (3.801d0*log10(T)) - (29487.d0/T))
  k11L = 10.d0**(-2.729d0 - (1.75d0*log10(T)) - (23474.d0/T))

  k8   = max(k8, 1d-40)
  k9   = max(k9, 1d-40)
  k10  = max(k10, 1d-40)
  k11  = max(k11, 1d-40)
  k9L  = max(k9L, 1d-40)
  k10L = max(k10L, 1d-40)
  k11L = max(k11L, 1d-40)

  !Log of all the rates
  lk8 = log10(k8)
  lk9 = (LTEfac*log10(k9L)) + (NLTEfac*log10(k9))
  lk10 = (LTEfac*log10(k10L)) + (NLTEfac*log10(k10))
  lk11 = (LTEfac*log10(k11L)) + (NLTEfac*log10(k11))

  rate = (ne*(10.d0**lk8)) + (nHI*(10.d0**lk9)) + (nH2*(10.d0**lk10)) + (nHeI*(10.0**lk11))
  rate = max(rate, 1d-40)

  rate = MAX(rate,1.d-100)

END FUNCTION beta_H2

FUNCTION alpha_CO(G0, xi_cr_H2, nCII, nH2, nOI) result(rate)
  ! see glover 2012
  implicit none

  real(dp), intent(in):: G0, xi_cr_H2, nCII, nH2, nOI
  real(dp):: rate
  real(dp):: k0, k1, gammaCHx_cr, gammaCHx, beta

  k0 = 5.d-16 ! cm^3 s^-1
  k1 = 5.d-10 ! Rate coefficient for the formation of CO from O + CHx

  ! Assuming CO formation is modulated by CH2+, https://home.strw.leidenuniv.nl/~ewine/photo/display_ch2+_65e31a07e69e64dbd64d37801983018f.html
  gammaCHx_cr = 8.88d-15 * (xi_cr_H2 / 1d-16) ! Cosmic rays
  gammaCHx = (1.41d-10 * G0 / G0_draine) + gammaCHx_cr ! 1.41e-10 per Draine field

  beta = k1 * nOI/(k1*nOI + gammaCHx)
  rate = k0 * nCII * nH2 * beta

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_CO

FUNCTION beta_CO(G0, xi_cr_H2) result(rate)
  ! CO destruction
  implicit none

  real(dp), intent(in):: G0, xi_cr_H2
  real(dp):: rate
  real(dp):: gammaCO, gammaCO_cr

  gammaCO = 2.43d-10 * G0 / G0_draine ! per Draine field (Heays+17)
  gammaCO_cr = 4.62d-15 * (xi_cr_H2 / 1d-16)

  rate = gammaCO + gammaCO_cr

  rate = MAX(rate,1.d-100)

END FUNCTION beta_CO

FUNCTION tau_dust_lw_mw(nH_tot, dx_SS, Z) result(tau)
  ! Dust optical depth in the Lyman-Werner band over dx_SS for MW dust scaled by Z (the dust-to-gas
  ! ratio in MW units), on the total hydrogen column (Draine & Bertoldi 1996)
  ! nH_tot : hydrogen nuclei number density, all states [cm^-3]
  implicit none
  real(dp), intent(in) :: nH_tot, dx_SS, Z
  real(dp) :: tau
  tau = sigma_dust_lw_mw * Z * nH_tot * dx_SS
END FUNCTION tau_dust_lw_mw

FUNCTION lw_transmission(tau, f_ani) result(trans)
  ! Mean-intensity transmission of a slab of normal optical depth tau for a field that is a beam
  ! (fraction f_ani, at normal incidence) plus an isotropic part: f e^-tau + (1-f) E2(tau). Exact in
  ! both limits, where exp(-(2-f) tau) only holds for tau << 1 (it is 4x too small at tau = 3).
  ! For 1 < tau <= 700 both terms take the same exp(-tau) (safe_exp(-tau) is that exp there).
  implicit none
  real(dp), intent(in) :: tau, f_ani
  real(dp) :: trans, f, ex
  f = min(max(f_ani, 0d0), 1d0)
  if (tau > 1d0 .and. tau <= 700d0) then
     ex = exp(-tau)
     trans = f * ex + (1d0 - f) * (ex * e2_scaled(tau))
  else
     trans = f * safe_exp(-tau) + (1d0 - f) * expint_E2(tau)
  end if
END FUNCTION lw_transmission

FUNCTION e2_scaled(x) result(g)
  ! e^x E2(x) for 1 <= x < 1024: the Chebyshev series of the sub-interval of x (E2_CHEB), by Clenshaw.
  ! x = y 2^e with y in [1, 2) (exponent, fraction: exact); sub-interval s = int(4 (y - 1)), local
  ! variable u = 8 (y - 1) - (2 s + 1) in [-1, 1) (exact)
  implicit none
  real(dp), intent(in) :: x
  real(dp) :: g, y, u, b0, b1, b2
  integer :: e, s, k, j
  e = exponent(x) - 1
  y = 2d0 * fraction(x)
  s = min(int(4d0 * (y - 1d0)), 3)
  k = min(max(4 * e + s, 0), E2_NI - 1)
  u = 8d0 * (y - 1d0) - dble(2 * s + 1)
  b1 = 0d0
  b2 = 0d0
  do j = E2_NC - 1, 1, -1
     b0 = 2d0 * u * b1 - b2 + E2_CHEB(j, k)
     b2 = b1
     b1 = b0
  end do
  g = u * b1 - b2 + E2_CHEB(0, k)
END FUNCTION e2_scaled

FUNCTION expint_E2(x) result(e2)
  ! Exponential integral E2(x), x >= 0: exp(-x) - x E1(x) with E1 by its series for x <= 1 (Numerical
  ! Recipes, expint); for x > 1, exp(-x) e2_scaled(x) (Chebyshev; relative error < 4e-16)
  implicit none
  real(dp), intent(in) :: x
  real(dp) :: e2, e1, term, fact
  real(dp), parameter :: euler = 0.5772156649015329d0, eps = 1d-14
  integer :: k
  if (x <= 0d0) then
     e2 = 1d0
     return
  end if
  if (x > 700d0) then
     e2 = 0d0
     return
  end if
  if (x <= 1d0) then
     e1 = -euler - log(x)
     fact = 1d0
     do k = 1, 100
        fact = -fact * x / k
        term = -fact / k
        e1 = e1 + term
        if (abs(term) < eps * abs(e1)) exit
     end do
  else
     e2 = exp(-x) * e2_scaled(x)
     return
  end if
  e2 = exp(-x) - x * e1
END FUNCTION expint_E2

FUNCTION comp_Sd(nH_tot, dx_SS, Z) result(ss_factor)
  ! Dust shielding of the Lyman-Werner band for MW dust scaled by Z, isotropic field
  ! (lw_transmission with f = 0); with CALIMA the per-bin grain cross-sections are used instead
  ! (compute_lw_dust_optical_depth in dust_photophysics.f90)
  implicit none
  real(dp), intent(in) :: nH_tot, dx_SS, Z
  real(dp) :: ss_factor
  ss_factor = lw_transmission(tau_dust_lw_mw(nH_tot, dx_SS, Z), 0d0)
END FUNCTION comp_Sd

FUNCTION comp_SH2(nH2, dx_SS) result(ss_factor)
  ! Returns the self shielding factor for dust
  ! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
  implicit none
  real(dp), intent(in):: nH2, dx_SS
  real(dp):: ss_factor
  real(dp):: xfac, cNH2, wH2, Sa, Sb, Sc

  cNH2 = nH2*dx_SS  !H2 column density
  xfac = cNH2/(5.d14)
  wH2 = 0.2d0

  Sa = (1.d0 - wH2)/((1.d0 + xfac)*(1.d0 + xfac))
  Sb = wH2/sqrt(1.d0 + xfac)
  Sc = safe_exp(-0.00085d0*sqrt(1.d0 + xfac))

  ss_factor = Sa + (Sb*Sc)
END FUNCTION comp_SH2

SUBROUTINE initialize_SCO_table()
  implicit none

  sco_table(:,1) = (/ &
  1.000d+00, 1.000d+12, 1.650d+12, 2.995d+12, 5.979d+12, 1.313d+13, &
  3.172d+13, 8.429d+13, 2.464d+14, 7.923d+14, 1.670d+15, 2.595d+15, &
  4.435d+15, 6.008d+15, 8.952d+15, 1.334d+16, 1.661d+16, 2.274d+16, &
  3.115d+16, 4.266d+16, 5.843d+16, 8.002d+16, 1.096d+17, 1.501d+17, &
  2.055d+17, 2.815d+17, 4.241d+17, 6.389d+17, 9.625d+17, 1.450d+18, &
  2.184d+18, 3.291d+18, 4.124d+18, 5.685d+18, 7.838d+18, 1.080d+19, &
  1.285d+19, 1.681d+19, 2.199d+19, 2.538d+19, 3.222d+19, 4.091d+19, &
  5.193d+19, 5.893d+19, 7.356d+19, 8.269d+19, 9.246d+19, 1.031d+20, &
  1.148d+20, 1.277d+20, 1.419d+20, 1.578d+20 /)

  sco_table(:,2) = (/ &
  1.000d+00, 9.990d-01, 9.981d-01, 9.961d-01, 9.912d-01, 9.815d-01, &
  9.601d-01, 9.113d-01, 8.094d-01, 6.284d-01, 4.808d-01, 3.889d-01, &
  2.827d-01, 2.293d-01, 1.695d-01, 1.224d-01, 1.017d-01, 7.764d-02, &
  5.931d-02, 4.546d-02, 3.506d-02, 2.728d-02, 2.143d-02, 1.700d-02, &
  1.360d-02, 1.094d-02, 8.273d-03, 6.283d-03, 4.773d-03, 3.611d-03, &
  2.704d-03, 1.986d-03, 1.657d-03, 1.258d-03, 9.332d-04, 6.745d-04, &
  5.596d-04, 4.123d-04, 2.982d-04, 2.490d-04, 1.827d-04, 1.324d-04, &
  9.473d-05, 7.891d-05, 5.668d-05, 4.732d-05, 3.967d-05, 3.327d-05, &
  2.788d-05, 2.331d-05, 1.944d-05, 1.619d-05 /)

  sh2_table(:,1) = (/ &
  1.000d+00, 2.666d+13, 3.801d+14, 6.634d+15, 8.829d+16, 9.268d+17, &
  1.007d+18, 2.021d+18, 3.036d+18, 4.051d+18, 5.066d+18, 6.082d+18, &
  7.097d+18, 8.112d+18, 9.341d+18, 1.014d+19, 2.030d+19, 3.045d+19, &
  4.061d+19, 5.076d+19, 6.092d+19, 7.107d+19, 8.123d+19, 9.353d+19, &
  1.015d+20, 2.031d+20, 3.047d+20, 4.062d+20, 5.078d+20, 6.094d+20, &
  7.109d+20, 8.125d+20, 9.355d+20, 1.016d+21, 2.031d+21, 3.047d+21, &
  4.063d+21, 5.078d+21, 6.094d+21, 7.110d+21, 8.125d+21, 9.355d+21, &
  1.016d+22 /) 

  sh2_table(:,2) = (/ &
  1.000d+00, 9.999d-01, 9.893d-01, 9.678d-01, 9.465d-01, 9.137d-01, &
  9.121d-01, 8.966d-01, 8.862d-01, 8.781d-01, 8.716d-01, 8.660d-01, &
  8.612d-01, 8.569d-01, 8.524d-01, 8.497d-01, 8.262d-01, 8.118d-01, &
  8.011d-01, 7.921d-01, 7.841d-01, 7.769d-01, 7.702d-01, 7.626d-01, &
  7.579d-01, 7.094d-01, 6.712d-01, 6.378d-01, 6.074d-01, 5.791d-01, &
  5.524d-01, 5.271d-01, 4.977d-01, 4.793d-01, 2.837d-01, 1.526d-01, &
  7.774d-02, 3.952d-02, 2.093d-02, 1.199d-02, 7.666d-03, 5.333d-03, &
  4.666d-03 /) 

END SUBROUTINE initialize_SCO_table

FUNCTION comp_SCO(nco_mol, nh2, dx_SS) result(ss_factor)
  ! Returns the self shielding factor for CO
  implicit none
  real(dp), intent(in):: nco_mol, nh2, dx_SS
  real(dp):: ss_factor
  real(dp):: logsCO, effcNCO, logeffcNCO
  real(dp):: logsH2, effcNH2, logeffcNH2
  integer:: i, idxCO, idxH2 ! Lower closest index

  ! initialize to 1.0
  ss_factor = 1.d0

  ! Pinned to table at lower boundary, extrapolated above upper boundary
  effcNCO = MAX(nco_mol*dx_SS, sco_table(1,1))
  if (effcNCO.ge.1.578d+20) then 
     ss_factor = ss_factor * 1.619d-05
  else
     ! Find the lower closest index
     idxCO = 1
     do i=1, 51 ! If above upper boundary then extrapolate using the slope 
        !between the last two points, so conveniently we set i to not go to 52
        if (effcNCO .ge. sco_table(i,1)) idxCO = i
     end do

     logeffcNCO = log10(effcNCO)
     ! Interpolate or extrapolate automatically
     logsCO = log10(sco_table(idxCO,2)) &
     + (log10(sco_table(idxCO+1,2)) - log10(sco_table(idxCO,2))) &
     * (logeffcNCO - log10(sco_table(idxCO,1))) &
     / (log10(sco_table(idxCO+1,1)) - log10(sco_table(idxCO,1)))

     ss_factor = ss_factor * (10.d0 ** logsCO)
  end if


  ! Pinned to table at lower boundary, extrapolated above upper boundary
  effcNH2 = MAX(nh2*dx_SS, sh2_table(1,1))
  if (effcNH2.ge.1.016d+22) then 
      ss_factor = ss_factor * 4.666d-03
  else 
      ! Find the lower closest index
      idxH2 = 1
      do i=1, 42 ! If above upper boundary then extrapolate using the slope 
        !between the last two points, so conveniently we set i to not go to 43
        if (effcNH2 .ge. sh2_table(i,1)) idxH2 = i
      end do

      logeffcNH2 = log10(effcNH2)
      ! Interpolate or extrapolate automatically
      logsH2 = log10(sh2_table(idxH2,2)) &
      + (log10(sh2_table(idxH2+1,2)) - log10(sh2_table(idxH2,2))) &
      * (logeffcNH2 - log10(sh2_table(idxH2,1))) &
      / (log10(sh2_table(idxH2+1,1)) - log10(sh2_table(idxH2,1)))

      ss_factor = ss_factor * (10.d0 ** logsH2)
  end if

END FUNCTION comp_SCO


end module molecules_module
