R1 = (ZZ/10007)[a, b, c, d, e, f, g]
I1 = ideal(
    f^3+4295*a^2*g-1316*a*b*g-3657*b^2*g-4810*a*c*g-4748*b*c*g+436*c^2*g-2744*a*d*g+1643*b*d*g+1794*c*d*g+2455*d^2*g+615*a*e*g-1008*b*e*g+154*c*e*g-3563*d*e*g-773*e^2*g+693*a*f*g+2331*b*f*g-4934*c*f*g-1817*d*f*g-1193*e*f*g-3493*f^2*g+3135*a*g^2+4693*b*g^2+2960*c*g^2+4957*d*g^2-392*e*g^2+874*f*g^2-4097*g^3,
    e*f^2+1144*a^2*g-2628*a*b*g+675*b^2*g+3183*a*c*g+1537*b*c*g-771*c^2*g-4529*a*d*g+4459*b*d*g-3927*c*d*g-936*d^2*g-3742*a*e*g-2110*b*e*g+4938*c*e*g+1395*d*e*g+527*e^2*g-3539*a*f*g+3375*b*f*g-3650*c*f*g-4522*d*f*g-2462*e*f*g-1146*f^2*g+4922*a*g^2+2571*b*g^2-417*c*g^2-217*d*g^2+2896*e*g^2+4368*f*g^2+3680*g^3,
    d*f^2+4484*a^2*g+2073*a*b*g+3134*b^2*g+3982*a*c*g-4054*b*c*g-535*c^2*g+3136*a*d*g+3820*b*d*g+2280*c*d*g+2262*d^2*g-4558*a*e*g+1419*b*e*g-4800*c*e*g-3496*d*e*g+4251*e^2*g+4187*a*f*g-492*b*f*g-2736*c*f*g+4798*d*f*g+3374*e*f*g-938*f^2*g+2768*a*g^2+1244*b*g^2+2645*c*g^2+1072*d*g^2+3694*e*g^2+2346*f*g^2-3181*g^3,
    c*f^2-2966*a^2*g+3904*a*b*g-4403*b^2*g+3143*a*c*g+1368*b*c*g+4230*c^2*g+1322*a*d*g-3366*b*d*g+217*c*d*g-406*d^2*g+989*a*e*g+3159*b*e*g+3266*c*e*g-885*d*e*g-462*e^2*g-3420*a*f*g-1616*b*f*g-2197*c*f*g+2751*d*f*g-839*e*f*g+3977*f^2*g+1489*a*g^2-4030*b*g^2-3309*c*g^2+3402*d*g^2-1971*e*g^2-3862*f*g^2+4887*g^3,
    b*f^2+2532*a^2*g+253*a*b*g+1983*b^2*g+3837*a*c*g-1887*b*c*g+3362*c^2*g-717*a*d*g+2523*b*d*g-3807*c*d*g-2534*d^2*g+556*a*e*g-2720*b*e*g-772*c*e*g+254*d*e*g+997*e^2*g+2163*a*f*g+1145*b*f*g+4104*c*f*g-1065*d*f*g-1719*e*f*g-2906*f^2*g-2990*a*g^2+4897*b*g^2-1957*c*g^2+3945*d*g^2-402*e*g^2+1311*f*g^2-4731*g^3,
    a*f^2+1410*a^2*g-1469*a*b*g-2214*b^2*g-638*a*c*g-256*b*c*g-3916*c^2*g+3895*a*d*g+2017*b*d*g+2009*c*d*g-2956*d^2*g-1210*a*e*g+3767*b*e*g+2035*c*e*g+591*d*e*g-781*e^2*g+4818*a*f*g-24*b*f*g+2155*c*f*g+4222*d*f*g+3366*e*f*g+3187*f^2*g+1250*a*g^2-2420*b*g^2-4304*c*g^2-1648*d*g^2-1413*e*g^2-147*f*g^2+1769*g^3,
    e^2*f-1687*a^2*g-743*a*b*g+2983*b^2*g-373*a*c*g+1781*b*c*g-4287*c^2*g-1864*a*d*g+2225*b*d*g-77*c*d*g-3538*d^2*g-3092*a*e*g+1949*b*e*g+4034*c*e*g-578*d*e*g+566*e^2*g-4093*a*f*g-84*b*f*g+354*c*f*g-1270*d*f*g+4491*e*f*g-2948*f^2*g-1362*a*g^2-2956*b*g^2-1762*c*g^2+1613*d*g^2-300*e*g^2-2519*f*g^2+197*g^3,
    d*e*f-4308*a^2*g+1652*a*b*g+4933*b^2*g-4361*a*c*g-4171*b*c*g-3839*c^2*g-174*a*d*g-2224*b*d*g+1162*c*d*g+2581*d^2*g+300*a*e*g-3150*b*e*g+1293*c*e*g-732*d*e*g+2580*e^2*g+2207*a*f*g+896*b*f*g-2533*c*f*g-518*d*f*g+4257*e*f*g-321*f^2*g+4771*a*g^2+4348*b*g^2+853*c*g^2+869*d*g^2-88*e*g^2-659*f*g^2-68*g^3,
    c*e*f+853*a^2*g-3946*a*b*g+1000*b^2*g-3453*a*c*g+2887*b*c*g-3673*c^2*g+1092*a*d*g+3832*b*d*g-1578*c*d*g-208*d^2*g+3207*a*e*g+1340*b*e*g+3639*c*e*g+1221*d*e*g-2285*e^2*g-797*a*f*g+344*b*f*g-1673*c*f*g-4153*d*f*g+2220*e*f*g-2907*f^2*g-452*a*g^2+2174*b*g^2-710*c*g^2-78*d*g^2-2829*e*g^2-2519*f*g^2-4325*g^3,
    b*e*f+37*a^2*g+801*a*b*g-4471*b^2*g-4396*a*c*g-4737*b*c*g-2919*c^2*g+4297*a*d*g+2539*b*d*g+4605*c*d*g+3669*d^2*g-4596*a*e*g-897*b*e*g+3598*c*e*g-1194*d*e*g+680*e^2*g+3903*a*f*g+1672*b*f*g-2158*c*f*g+119*d*f*g-4215*e*f*g-3513*f^2*g+1978*a*g^2+4350*b*g^2-3629*c*g^2-4444*d*g^2-2962*e*g^2-3950*f*g^2-2643*g^3,
    a*e*f-1495*a^2*g-1380*a*b*g+4660*b^2*g+404*a*c*g+1956*b*c*g-221*c^2*g-3844*a*d*g-3047*b*d*g-2386*c*d*g+1994*d^2*g-3621*a*e*g+3131*b*e*g+2834*c*e*g+4836*d*e*g+577*e^2*g+4499*a*f*g+592*b*f*g-4847*c*f*g+1551*d*f*g-1525*e*f*g-3269*f^2*g-4455*a*g^2-1534*b*g^2-1006*c*g^2+3681*d*g^2+4780*e*g^2+4922*f*g^2+4236*g^3,
    d^2*f-4939*a^2*g-1847*a*b*g+4618*b^2*g+2137*a*c*g+3327*b*c*g+3008*c^2*g-1888*a*d*g+1093*b*d*g+88*c*d*g+3061*d^2*g+4997*a*e*g-3152*b*e*g-1387*c*e*g+3535*d*e*g+666*e^2*g+1356*a*f*g-3221*b*f*g+3444*c*f*g-395*d*f*g-4038*e*f*g+754*f^2*g-1880*a*g^2+3192*b*g^2-3295*c*g^2-256*d*g^2+988*e*g^2+3745*f*g^2-4176*g^3,
    c*d*f-175*a^2*g-3247*a*b*g-1681*b^2*g-4659*a*c*g+1997*b*c*g-511*c^2*g+3186*a*d*g-2010*b*d*g+741*c*d*g-1053*d^2*g+3814*a*e*g+1673*b*e*g-607*c*e*g+4205*d*e*g-4044*e^2*g+43*a*f*g-4094*b*f*g-517*c*f*g-2407*d*f*g-3202*e*f*g+1110*f^2*g-2974*a*g^2+4997*b*g^2-4014*c*g^2+4583*d*g^2+1265*e*g^2+2821*f*g^2-1406*g^3,
    b*d*f+890*a^2*g+430*a*b*g-3688*b^2*g-3783*a*c*g+4034*b*c*g-4901*c^2*g-3338*a*d*g+3945*b*d*g+713*c*d*g-470*d^2*g-2791*a*e*g+3276*b*e*g-1262*c*e*g+2572*d*e*g+226*e^2*g+1572*a*f*g-3919*b*f*g+1849*c*f*g-1845*d*f*g-4383*e*f*g+1627*f^2*g+338*a*g^2+2303*b*g^2+3373*c*g^2-3606*d*g^2+1257*e*g^2+1288*f*g^2+3800*g^3,
    a*d*f+3862*a^2*g-3638*a*b*g-3612*b^2*g+4073*a*c*g-4675*b*c*g+3885*c^2*g+1051*a*d*g+3196*b*d*g+2172*c*d*g-322*d^2*g-111*a*e*g+458*b*e*g+2777*c*e*g+4235*d*e*g+4397*e^2*g+186*a*f*g-3664*b*f*g-2425*c*f*g-4730*d*f*g-4556*e*f*g-4525*f^2*g+4577*a*g^2+2044*b*g^2+2762*c*g^2-454*d*g^2-483*e*g^2-1850*f*g^2+1655*g^3,
    c^2*f+247*a^2*g-1644*a*b*g-2830*b^2*g-2308*a*c*g-3616*b*c*g-2822*c^2*g+2002*a*d*g+4149*b*d*g+2787*c*d*g+2205*d^2*g+2428*a*e*g+3085*b*e*g+141*c*e*g-2364*d*e*g+3247*e^2*g-598*a*f*g-389*b*f*g-1069*c*f*g-483*d*f*g-3829*e*f*g+679*f^2*g+1822*a*g^2+4704*b*g^2-3397*c*g^2+893*d*g^2-1800*e*g^2-4042*f*g^2-2119*g^3,
    b*c*f-4414*a^2*g-506*a*b*g-4944*b^2*g+446*a*c*g-4084*b*c*g+3920*c^2*g-1601*a*d*g+2458*b*d*g-3094*c*d*g+3089*d^2*g-836*a*e*g-133*b*e*g+3534*c*e*g-2565*d*e*g+2163*e^2*g+4258*a*f*g-644*b*f*g+1237*c*f*g+1467*d*f*g-3403*e*f*g-3357*f^2*g+2383*a*g^2-4573*b*g^2-4955*c*g^2+3585*d*g^2-704*e*g^2+4422*f*g^2+4692*g^3,
    a*c*f+2120*a^2*g-4271*a*b*g+845*b^2*g-1800*a*c*g+1562*b*c*g+793*c^2*g+2303*a*d*g-3641*b*d*g-951*c*d*g-1992*d^2*g+4534*a*e*g-3188*b*e*g-4338*c*e*g-4738*d*e*g-4467*e^2*g+892*a*f*g+1501*b*f*g+4288*c*f*g-1125*d*f*g-3442*e*f*g+4601*f^2*g-13*a*g^2-4788*b*g^2+1468*c*g^2-3360*d*g^2-4656*e*g^2-4570*f*g^2+1444*g^3,
    b^2*f+4911*a^2*g+4538*a*b*g+4709*b^2*g+1945*a*c*g-3666*b*c*g-1512*c^2*g+2964*a*d*g+3611*b*d*g-4531*c*d*g+347*d^2*g-4121*a*e*g-340*b*e*g-786*c*e*g+1501*d*e*g-1363*e^2*g+1044*a*f*g+2780*b*f*g+827*c*f*g+370*d*f*g+3994*e*f*g+2266*f^2*g+14*a*g^2-361*b*g^2+1809*c*g^2+847*d*g^2+1699*e*g^2+2842*f*g^2-4417*g^3,
    a*b*f-3682*a^2*g+2696*a*b*g-3057*b^2*g+2737*a*c*g-308*b*c*g-4226*c^2*g-3067*a*d*g-2110*b*d*g-4702*c*d*g-2440*d^2*g+1106*a*e*g-2814*b*e*g+2723*c*e*g+4549*d*e*g+2333*e^2*g-2097*a*f*g+3192*b*f*g-4012*c*f*g-4683*d*f*g+965*e*f*g-2262*f^2*g-4802*a*g^2-885*b*g^2-4269*c*g^2-2251*d*g^2+922*e*g^2-3028*f*g^2-4475*g^3,
    a^2*f-1430*a^2*g-4419*a*b*g+2829*b^2*g-4854*a*c*g+1760*b*c*g-2910*c^2*g+3194*a*d*g+4948*b*d*g-60*c*d*g-1697*d^2*g+3140*a*e*g-4183*b*e*g+4209*c*e*g+4055*d*e*g+3031*e^2*g+1197*a*f*g+3458*b*f*g-2045*c*f*g+4078*d*f*g+1262*e*f*g-3773*f^2*g+1284*a*g^2+4994*b*g^2-976*c*g^2+3357*d*g^2+4502*e*g^2-148*f*g^2+3141*g^3,
    e^3+3195*a^2*g+4072*a*b*g+2623*b^2*g+3649*a*c*g-3372*b*c*g+1415*c^2*g-2533*a*d*g-4870*b*d*g-1969*c*d*g+948*d^2*g-4495*a*e*g-3236*b*e*g+4601*c*e*g+890*d*e*g-3292*e^2*g+899*a*f*g+798*b*f*g-512*c*f*g-1928*d*f*g+201*e*f*g-751*f^2*g-3239*a*g^2+4087*b*g^2+4769*c*g^2+2332*d*g^2+652*e*g^2+4645*f*g^2+942*g^3,
    d*e^2+4751*a^2*g+1439*a*b*g+3943*b^2*g-3574*a*c*g-915*b*c*g+2908*c^2*g-4775*a*d*g-2890*b*d*g-4076*c*d*g-4703*d^2*g+961*a*e*g-1850*b*e*g+4632*c*e*g-2011*d*e*g-3633*e^2*g-1709*a*f*g-3581*b*f*g+4510*c*f*g-4662*d*f*g-2571*e*f*g-1795*f^2*g+1288*a*g^2-2873*b*g^2-973*c*g^2+2826*d*g^2-3222*e*g^2+1562*f*g^2-3816*g^3,
    c*e^2+2603*a^2*g-1277*a*b*g+1776*b^2*g+1701*a*c*g+2528*b*c*g+2150*c^2*g-1880*a*d*g-740*b*d*g+3968*c*d*g+3503*d^2*g+1814*a*e*g-4327*b*e*g-2719*c*e*g+1029*d*e*g+4282*e^2*g-1186*a*f*g+4740*b*f*g-521*c*f*g+2554*d*f*g-4387*e*f*g-4469*f^2*g-4504*a*g^2-111*b*g^2-363*c*g^2+3327*d*g^2+2091*e*g^2+2206*f*g^2+4880*g^3,
    b*e^2-2100*a^2*g-4403*a*b*g+3223*b^2*g-407*a*c*g-1005*b*c*g+441*c^2*g+1154*a*d*g-765*b*d*g+3249*c*d*g+4773*d^2*g+4671*a*e*g+785*b*e*g+4914*c*e*g+983*d*e*g-100*e^2*g+3544*a*f*g+585*b*f*g+636*c*f*g+1358*d*f*g+1336*e*f*g-145*f^2*g+1172*a*g^2-665*b*g^2+2339*c*g^2-4336*d*g^2+1163*e*g^2+1658*f*g^2-4979*g^3,
    a*e^2-2696*a^2*g-916*a*b*g+903*b^2*g+2699*a*c*g-4676*b*c*g-393*c^2*g+2464*a*d*g+1044*b*d*g-2382*c*d*g-1735*d^2*g+929*a*e*g-4886*b*e*g+4312*c*e*g+3165*d*e*g-1976*e^2*g+1681*a*f*g+3061*b*f*g+822*c*f*g-4400*d*f*g-2307*e*f*g-697*f^2*g+4372*a*g^2-383*b*g^2+104*c*g^2+1343*d*g^2+4921*e*g^2+1541*f*g^2+3042*g^3,
    d^2*e+2515*a^2*g+2141*a*b*g+1263*b^2*g+515*a*c*g+3584*b*c*g+2611*c^2*g+816*a*d*g+465*b*d*g+2234*c*d*g-3987*d^2*g-3619*a*e*g-3611*b*e*g+1245*c*e*g+2748*d*e*g-1128*e^2*g-a*f*g+1607*b*f*g+926*c*f*g+2130*d*f*g-2306*e*f*g-1965*f^2*g+2090*a*g^2+3669*b*g^2+3120*c*g^2+3273*d*g^2+4168*e*g^2+3821*f*g^2-375*g^3,
    c*d*e+4541*a^2*g+891*a*b*g-3742*b^2*g-3931*a*c*g+861*b*c*g+4920*c^2*g+2687*a*d*g+4590*b*d*g-497*c*d*g+1728*d^2*g+1266*a*e*g+2830*b*e*g-2179*c*e*g+3399*d*e*g-1284*e^2*g+1500*a*f*g-4954*b*f*g+2367*c*f*g-155*d*f*g+2253*e*f*g+3650*f^2*g-3736*a*g^2-2586*b*g^2+1793*c*g^2+4976*d*g^2+2407*e*g^2-647*f*g^2-1995*g^3,
    b*d*e-4707*a^2*g+2612*a*b*g+1259*b^2*g-1052*a*c*g-2301*b*c*g+3591*c^2*g-3620*a*d*g-1624*b*d*g+2140*c*d*g-476*d^2*g+2925*a*e*g+3898*b*e*g+380*c*e*g-66*d*e*g-708*e^2*g+3920*a*f*g+205*b*f*g+15*c*f*g+3675*d*f*g-4015*e*f*g-4275*f^2*g-3483*a*g^2-3076*b*g^2+2893*c*g^2-19*d*g^2-4327*e*g^2-46*f*g^2+2573*g^3,
    a*d*e-13*a^2*g+2340*a*b*g-210*b^2*g-3963*a*c*g+3029*b*c*g+4409*c^2*g+3984*a*d*g+2501*b*d*g+432*c*d*g-2465*d^2*g+3424*a*e*g+1939*b*e*g-2349*c*e*g+4776*d*e*g+2297*e^2*g-1868*a*f*g+1126*b*f*g+4209*c*f*g-731*d*f*g-3253*e*f*g+2506*f^2*g-3132*a*g^2-4101*b*g^2+22*c*g^2+67*d*g^2+1718*e*g^2+1876*f*g^2+3889*g^3,
    c^2*e+2845*a^2*g+2879*a*b*g+2178*b^2*g-4703*a*c*g-1727*b*c*g+3912*c^2*g+457*a*d*g-1711*b*d*g-869*c*d*g-2791*d^2*g+1582*a*e*g+1666*b*e*g-717*c*e*g+4876*d*e*g-4533*e^2*g-789*a*f*g+1949*b*f*g-2357*c*f*g+3517*d*f*g+1441*e*f*g-4647*f^2*g-2119*a*g^2+1573*b*g^2+777*c*g^2-586*d*g^2-3094*e*g^2+3265*f*g^2-1541*g^3,
    b*c*e-4615*a^2*g-3734*a*b*g+863*b^2*g-2317*a*c*g-3589*b*c*g+2720*c^2*g-2941*a*d*g-2941*b*d*g-1223*c*d*g-4108*d^2*g+3372*a*e*g+345*b*e*g-1206*c*e*g+2314*d*e*g+3382*e^2*g-1646*a*f*g-1286*b*f*g+845*c*f*g+8*d*f*g-3268*e*f*g+1316*f^2*g-1033*a*g^2-881*b*g^2+1302*c*g^2-4254*d*g^2+3852*e*g^2-3643*f*g^2+987*g^3,
    a*c*e-1834*a^2*g+3684*a*b*g+3952*b^2*g-797*a*c*g+1549*b*c*g-1010*c^2*g+4277*a*d*g+1546*b*d*g+1371*c*d*g-2394*d^2*g+2823*a*e*g-2474*b*e*g+2914*c*e*g-1490*d*e*g+1220*e^2*g-4572*a*f*g+77*b*f*g+180*c*f*g+694*d*f*g-2832*e*f*g-3036*f^2*g+2495*a*g^2-4638*b*g^2-4504*c*g^2+1239*d*g^2+2518*e*g^2+2445*f*g^2+3284*g^3,
    b^2*e-2393*a^2*g-4190*a*b*g-3900*b^2*g+1594*a*c*g-1119*b*c*g-2041*c^2*g-1825*a*d*g+1838*b*d*g-4965*c*d*g+4281*d^2*g-2453*a*e*g+497*b*e*g+566*c*e*g+2068*d*e*g-3477*e^2*g+1759*a*f*g+1398*b*f*g-1451*c*f*g-43*d*f*g+4115*e*f*g-89*f^2*g+299*a*g^2+3148*b*g^2+269*c*g^2-2744*d*g^2+1159*e*g^2-1123*f*g^2+1615*g^3,
    a*b*e-772*a^2*g-4938*a*b*g+4036*b^2*g-2204*a*c*g+2841*b*c*g+3338*c^2*g-1428*a*d*g-2536*b*d*g+1826*c*d*g-2579*d^2*g+294*a*e*g-4412*b*e*g-2861*c*e*g+3074*d*e*g+1014*e^2*g-592*a*f*g-878*b*f*g-224*c*f*g+4531*d*f*g-4175*e*f*g+1815*f^2*g-4366*a*g^2-1826*b*g^2+2471*c*g^2+4267*d*g^2-4226*e*g^2-1081*f*g^2+2553*g^3,
    a^2*e+4703*a^2*g-3675*a*b*g-794*b^2*g+2851*a*c*g+2021*b*c*g-1696*c^2*g+3112*a*d*g-2774*b*d*g-4489*c*d*g+1726*d^2*g-4281*a*e*g-3831*b*e*g+3962*c*e*g-4437*d*e*g+4955*e^2*g-534*a*f*g-1038*b*f*g+4056*c*f*g+3272*d*f*g-4686*e*f*g+1017*f^2*g-1787*a*g^2+4250*b*g^2-4638*c*g^2+2514*d*g^2+4309*e*g^2+4757*f*g^2-4281*g^3,
    d^3+2719*a^2*g-2847*a*b*g+2407*b^2*g+4350*a*c*g+2049*b*c*g-1772*c^2*g-2046*a*d*g-3580*b*d*g-822*c*d*g+4435*d^2*g-2785*a*e*g+3960*b*e*g+1408*c*e*g+2768*d*e*g+3396*e^2*g-1034*a*f*g+3922*b*f*g+1759*c*f*g-4616*d*f*g-709*e*f*g+2679*f^2*g-2766*a*g^2+2847*b*g^2+1119*c*g^2+501*d*g^2+3339*e*g^2-1125*f*g^2-334*g^3,
    c*d^2-3038*a^2*g+123*a*b*g-2875*b^2*g-426*a*c*g+1502*b*c*g-1673*c^2*g+4233*a*d*g-273*b*d*g+4239*c*d*g-2462*d^2*g-4034*a*e*g-1188*b*e*g+3879*c*e*g-3263*d*e*g-3191*e^2*g+361*a*f*g+3450*b*f*g-1307*c*f*g-315*d*f*g-1315*e*f*g+4530*f^2*g+2061*a*g^2-2919*b*g^2-1471*c*g^2+3789*d*g^2-1860*e*g^2+1479*f*g^2+4813*g^3,
    b*d^2+3691*a^2*g+3400*a*b*g-61*b^2*g-587*a*c*g-2618*b*c*g+4483*c^2*g-1135*a*d*g-3399*b*d*g-3677*c*d*g+852*d^2*g+4395*a*e*g-614*b*e*g-4132*c*e*g+2357*d*e*g-1294*e^2*g+1987*a*f*g+3147*b*f*g-4748*c*f*g-2011*d*f*g+3479*e*f*g-510*f^2*g+1584*a*g^2+4231*b*g^2-241*c*g^2+879*d*g^2-2071*e*g^2+2870*f*g^2+3974*g^3,
    a*d^2-3421*a^2*g+3831*a*b*g+1839*b^2*g+3847*a*c*g-341*b*c*g-1981*c^2*g-3935*a*d*g-2160*b*d*g-4473*c*d*g-3902*d^2*g-3829*a*e*g-3413*b*e*g-2363*c*e*g-4146*d*e*g-481*e^2*g-2300*a*f*g+4245*b*f*g-1533*c*f*g-4327*d*f*g-1629*e*f*g-1792*f^2*g+2692*a*g^2+2185*b*g^2-4844*c*g^2-4110*d*g^2-1308*e*g^2+3342*f*g^2-4240*g^3,
    c^2*d-1995*a^2*g-4763*a*b*g-3616*b^2*g+445*a*c*g+1383*b*c*g+4078*c^2*g-2748*a*d*g-2562*b*d*g-2068*c*d*g-3006*d^2*g+1502*a*e*g-1713*b*e*g+1347*c*e*g+4739*d*e*g-1426*e^2*g+1543*a*f*g+1547*b*f*g+4925*c*f*g-3767*d*f*g-286*e*f*g+1081*f^2*g+3153*a*g^2-4321*b*g^2+1123*c*g^2-2613*d*g^2+3637*e*g^2+1994*f*g^2-4506*g^3,
    b*c*d+4575*a^2*g+3076*a*b*g+4026*b^2*g+4369*a*c*g+4654*b*c*g-2331*c^2*g+2784*a*d*g-2858*b*d*g-3884*c*d*g-2175*d^2*g-4386*a*e*g+4904*b*e*g-3393*c*e*g+3304*d*e*g+4190*e^2*g+3942*a*f*g+4869*b*f*g-4423*c*f*g-3781*d*f*g-1125*e*f*g-2730*f^2*g+3440*a*g^2-470*b*g^2-2928*c*g^2-328*d*g^2+2096*e*g^2+3543*f*g^2-1609*g^3,
    a*c*d-23*a^2*g+3*a*b*g-1614*b^2*g-4379*a*c*g+3020*b*c*g+4822*c^2*g+4644*a*d*g+4278*b*d*g-1817*c*d*g+1517*d^2*g-1175*a*e*g+4904*b*e*g+4410*c*e*g-841*d*e*g+1081*e^2*g-4406*a*f*g+2085*b*f*g+4200*c*f*g-1927*d*f*g-3437*e*f*g+4422*f^2*g-4789*a*g^2-3110*b*g^2+3194*c*g^2+4752*d*g^2+1713*e*g^2+2117*f*g^2+362*g^3,
    b^2*d-982*a^2*g+2915*a*b*g-2358*b^2*g-972*a*c*g-206*b*c*g+495*c^2*g+4854*a*d*g+1527*b*d*g-4280*c*d*g+2685*d^2*g-2505*a*e*g+1727*b*e*g-4986*c*e*g+1232*d*e*g-2178*e^2*g+1891*a*f*g-1596*b*f*g-3878*c*f*g-112*d*f*g-1363*e*f*g-2366*f^2*g+3301*a*g^2+3669*b*g^2+3838*c*g^2+3820*d*g^2-439*e*g^2-464*f*g^2+4809*g^3,
    a*b*d+4491*a^2*g+289*a*b*g+4009*b^2*g-4209*a*c*g+2389*b*c*g-1986*c^2*g-2138*a*d*g-3321*b*d*g+647*c*d*g+4221*d^2*g+3610*a*e*g-3394*b*e*g-2370*c*e*g+212*d*e*g+3011*e^2*g+1084*a*f*g+2448*b*f*g-3456*c*f*g+129*d*f*g-1561*e*f*g-4095*f^2*g-4732*a*g^2+4532*b*g^2+1560*c*g^2-4022*d*g^2+3294*e*g^2-4196*f*g^2+1890*g^3,
    a^2*d+2585*a^2*g+523*a*b*g-4286*b^2*g+3948*a*c*g-2577*b*c*g+1142*c^2*g+3001*a*d*g-4912*b*d*g-4300*c*d*g+605*d^2*g+2980*a*e*g+2309*b*e*g-2767*c*e*g-2515*d*e*g+2280*e^2*g-2015*a*f*g+1361*b*f*g+2237*c*f*g-383*d*f*g+2084*e*f*g-1416*f^2*g-3165*a*g^2-1891*b*g^2+3443*c*g^2-2659*d*g^2-1886*e*g^2+1309*f*g^2-2607*g^3,
    c^3+1848*a^2*g-1396*a*b*g+4221*b^2*g+1884*a*c*g-214*b*c*g-4022*c^2*g+557*a*d*g+1810*b*d*g+445*c*d*g-1178*d^2*g-1247*a*e*g-4035*b*e*g+3612*c*e*g+1879*d*e*g-3411*e^2*g-386*a*f*g+4680*b*f*g+4226*c*f*g+810*d*f*g-3283*e*f*g-3916*f^2*g+1062*a*g^2-645*b*g^2+3123*c*g^2+2743*d*g^2-3835*e*g^2-3295*f*g^2+1173*g^3,
    b*c^2-3704*a^2*g-2703*a*b*g+313*b^2*g-4745*a*c*g+4516*b*c*g+108*c^2*g-3867*a*d*g+1237*b*d*g+1996*c*d*g+155*d^2*g-2467*a*e*g-3500*b*e*g-2090*c*e*g+3684*d*e*g-3518*e^2*g-3720*a*f*g+2203*b*f*g-862*c*f*g+2168*d*f*g-3285*e*f*g-418*f^2*g+2882*a*g^2+1606*b*g^2+469*c*g^2+1605*d*g^2-4794*e*g^2+3866*f*g^2+75*g^3,
    a*c^2+2348*a^2*g-3936*a*b*g+4319*b^2*g+442*a*c*g+3244*b*c*g-3607*c^2*g-244*a*d*g-2126*b*d*g-2048*c*d*g-1137*d^2*g+91*a*e*g-1220*b*e*g-2381*c*e*g+1167*d*e*g-1306*e^2*g-400*a*f*g+3379*b*f*g+3654*c*f*g-3696*d*f*g-2800*e*f*g+1110*f^2*g-3595*a*g^2-3792*b*g^2-3971*c*g^2-180*d*g^2-1512*e*g^2-4740*f*g^2+3863*g^3,
    b^2*c+3212*a^2*g+1288*a*b*g-2523*b^2*g+4465*a*c*g+3636*b*c*g-3500*c^2*g-2024*a*d*g-2125*b*d*g-2537*c*d*g-2406*d^2*g+2890*a*e*g+3553*b*e*g-213*c*e*g+4957*d*e*g+859*e^2*g+2690*a*f*g+4983*b*f*g-3408*c*f*g+4418*d*f*g-315*e*f*g+711*f^2*g+751*a*g^2+1759*b*g^2-4949*c*g^2+4947*d*g^2-4484*e*g^2+1709*f*g^2-4135*g^3,
    a*b*c+2833*a^2*g+3169*a*b*g-2402*b^2*g-2364*a*c*g-990*b*c*g-1746*c^2*g+4901*a*d*g+4039*b*d*g+2872*c*d*g+146*d^2*g-3416*a*e*g+2361*b*e*g-143*c*e*g-1147*d*e*g+64*e^2*g+222*a*f*g-1126*b*f*g-4948*c*f*g-371*d*f*g-4059*e*f*g+2759*f^2*g-2769*a*g^2-2270*b*g^2-3765*c*g^2+4770*d*g^2+2974*e*g^2+3747*f*g^2+2798*g^3,
    a^2*c-4706*a^2*g-3060*a*b*g-4084*b^2*g+1381*a*c*g-2910*b*c*g-3170*c^2*g+3014*a*d*g+2639*b*d*g+3654*c*d*g+3854*d^2*g-2426*a*e*g+2920*b*e*g-2882*c*e*g-3953*d*e*g-4801*e^2*g-201*a*f*g+4617*b*f*g+3527*c*f*g-4060*d*f*g+2425*e*f*g+2538*f^2*g+214*a*g^2-2283*b*g^2+4812*c*g^2+288*d*g^2-3477*e*g^2-1647*f*g^2+4270*g^3,
    b^3+892*a^2*g-3732*a*b*g+1959*b^2*g-1569*a*c*g-1340*b*c*g+58*c^2*g+2391*a*d*g-909*b*d*g-4164*c*d*g-506*d^2*g-880*a*e*g+4757*b*e*g-174*c*e*g-4179*d*e*g+5000*e^2*g-405*a*f*g-2381*b*f*g+4855*c*f*g+2589*d*f*g-4236*e*f*g-3523*f^2*g-2165*a*g^2+2640*b*g^2+2615*c*g^2-3992*d*g^2+1749*e*g^2+2690*f*g^2+582*g^3,
    a*b^2-2369*a^2*g+1845*a*b*g+3164*b^2*g+1173*a*c*g-2801*b*c*g+4410*c^2*g+402*a*d*g+1532*b*d*g-2089*c*d*g-903*d^2*g+4803*a*e*g+810*b*e*g-557*c*e*g+805*d*e*g-4465*e^2*g-915*a*f*g+3963*b*f*g-3859*c*f*g-4516*d*f*g-1271*e*f*g-4979*f^2*g-3036*a*g^2+4867*b*g^2+4431*c*g^2-1848*d*g^2-2326*e*g^2+2029*f*g^2-4133*g^3,
    a^2*b-2348*a^2*g-586*a*b*g+1731*b^2*g-3150*a*c*g-3068*b*c*g-1014*c^2*g+817*a*d*g-2514*b*d*g-3811*c*d*g+4341*d^2*g-1827*a*e*g+1748*b*e*g-4987*c*e*g-2464*d*e*g-3219*e^2*g-3494*a*f*g+2211*b*f*g-3954*c*f*g+3753*d*f*g-2110*e*f*g-2034*f^2*g+3745*a*g^2-2231*b*g^2-3941*c*g^2+3847*d*g^2-4041*e*g^2-3452*f*g^2+630*g^3,
    a^3-2598*a^2*g+432*a*b*g-3501*b^2*g+3941*a*c*g+1831*b*c*g+2605*c^2*g+22*a*d*g-3375*b*d*g-1660*c*d*g+614*d^2*g-478*a*e*g-1745*b*e*g+784*c*e*g+2721*d*e*g+1032*e^2*g-1323*a*f*g-1900*b*f*g-927*c*f*g-1298*d*f*g-2745*e*f*g+387*f^2*g+2505*a*g^2+4836*b*g^2+300*c*g^2+2250*d*g^2+1067*e*g^2-4535*f*g^2-4784*g^3
    )

