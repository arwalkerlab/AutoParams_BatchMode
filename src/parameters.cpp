#include "parameters.h"


std::map <std::string, std::string> BASIC_ATOM_TYPES={
//sp carbons
{"CZ","C_sp"}, {"CY","C_sp"}, {"cg","C_sp"}, {"ch","C_sp"}, {"c1","C_sp"},
//sp2 carbons
{"C4","C_sp2"}, {"C5","C_sp2"}, {"C ","C_sp2"}, {"CA","C_sp2"}, {"CB","C_sp2"},
{"CC","C_sp2"}, {"CD","C_sp2"}, {"CK","C_sp2"}, {"CP","C_sp2"}, {"CM","C_sp2"},
{"CS","C_sp2"}, {"CN","C_sp2"}, {"CQ","C_sp2"}, {"CR","C_sp2"}, {"CV","C_sp2"},
{"CW","C_sp2"}, {"C*","C_sp2"}, {"CI","C_sp2"}, {"Ck","C_sp2"}, {"Cj","C_sp2"},
{"cC","C_sp2"}, {"cB","C_sp2"}, {"c ","C_sp2"}, {"ce","C_sp2"}, {"cf","C_sp2"},
{"ca","C_sp2"},
//sp3 carbons
{"CX","C_sp3"}, {"CT","C_sp3"}, {"cA","C_sp3"}, {"cD","C_sp3"}, {"cP","C_sp3"}, 
{"cR","C_sp3"}, {"Cg","C_sp3"}, {"Cy","C_sp3"}, {"2C","C_sp3"}, {"3C","C_sp3"}, 
{"Cp","C_sp3"}, {"cx","C_sp3"}, {"cy","C_sp3"}, {"CG","C_sp3"},
//sp nitrogens
{"NY","N_sp"}, {"n1","N_sp"}, {"no","N_sp"},
//sp2 nitrogens
{"N ","N_sp2"}, {"NA","N_sp2"}, {"NB","N_sp2"}, {"NC","N_sp2"}, {"N2","N_sp2"},
{"N*","N_sp2"}, {"nN","N_sp2"}, {"Ng","N_sp2"}, {"nb","N_sp2"}, {"nc","N_sp2"},
{"nd","N_sp2"}, {"n ","N_sp2"}, {"ni","N_sp2"}, {"nj","N_sp2"}, {"ne","N_sp2"},
{"nf","N_sp2"}, {"n2","N_sp2"}, 
//sp3 nitrogens
{"N3","N_sp3"}, {"NT","N_sp3"}, {"nA","N_sp3"}, {"n4","N_sp3"}, {"nk","N_sp3"},
{"nl","N_sp3"}, {"np","N_sp3"}, {"n3","N_sp3"}, {"nq","N_sp3"}, {"nh","N_sp3"},
{"nm","N_sp3"}, {"nn","N_sp3"}, {"na","N_sp3"}, 
//single-bond oxygens
{"OW","O_sb"}, {"OH","O_sb"}, {"OS","O_sb"}, {"oH","O_sb"}, {"oR","O_sb"},
{"oS","O_sb"}, {"oT","O_sb"}, {"Os","O_sb"}, {"Oy","O_sb"}, {"Oh","O_sb"},
{"ow","O_sb"}, {"oh","O_sb"}, {"os","O_sb"}, {"oq","O_sb"}, 
//double-bond oxygens
{"OP","O_db"}, {"O2","O_db"}, {"O ","O_db"}, {"oC","O_db"}, {"oO","O_db"},
{"oP","O_db"}, {"o ","O_db"}, 
//hydrogens
{"HW","H"}, {"H4","H"}, {"H5","H"}, {"HA","H"}, {"HC","H"}, {"H1","H"}, 
{"H2","H"}, {"H3","H"}, {"HP","H"}, {"HZ","H"}, {"H ","H"}, {"HO","H"},
{"HS","H"}, {"hA","H"}, {"hL","H"}, {"hE","H"}, {"hS","H"}, {"hX","H"},
{"hB","H"}, {"hN","H"}, {"hO","H"}, {"hR","H"}, {"Hc","H"}, {"Hp","H"},
{"Ha","H"}, {"Ho","H"}, {"hc","H"}, {"h1","H"}, {"h2","H"}, {"h3","H"},
{"hx","H"}, {"ha","H"}, {"h4","H"}, {"h5","H"}, {"hn","H"}, {"ho","H"},
{"hs","H"}, {"hw","H"}, {"hp","H"},
//halogens
{"F ","F"}, {"f ","F"},
{"Cl","Cl"}, {"cl","Cl"}, {"CL","Cl"},
{"Br","Br"}, {"br","Br"}, {"BR","Br"},
{"I ","I"}, {"i ","I"},
//sulfur
{"S ","S"}, {"SH","S"}, {"SO","S"}, {"SS","S"}, {"Sm","S"}, {"SF","S"},
{"S*","S"}, {"s ","S"}, {"s2","S"}, {"s4","S"}, {"s6","S"}, {"sh","S"},
{"ss","S"}, {"sp","S"}, {"sq","S"}, {"sx","S"}, {"sy","S"}, {"SG","S"},
//phosphorus
{"P ","P"}, {"pA","P"}, {"p2","P"}, {"p3","P"}, {"p4","P"}, {"p5","P"},
{"pb","P"}, {"pc","P"}, {"pd","P"}, {"pe","P"}, {"pf","P"}, {"px","P"},
{"py","P"},
//lone pair.
{"EP","EP"}};

// Masses library (map of bond-pair strings to full frcmod string)
std::map<std::string,std::string> mass_library={
{"C ","C  12.01         0.616"},{"CA","CA 12.01         0.360"},{"CB","CB 12.01         0.360"},{"CC","CC 12.01         0.360"},{"CD","CD 12.01         0.360"},
{"CI","CI 12.01              "},{"CK","CK 12.01         0.360"},{"CP","CP 12.01         0.360"},{"CM","CM 12.01         0.360"},{"CS","CS 12.01         0.360"},
{"CN","CN 12.01         0.360"},{"CQ","CQ 12.01         0.360"},{"CR","CR 12.01         0.360"},{"CT","CT 12.01         0.878"},{"CV","CV 12.01         0.360"},
{"CW","CW 12.01         0.360"},{"C*","C* 12.01         0.360"},{"CX","CX 12.01         0.360"},{"CY","CY 12.01         0.360"},{"CZ","CZ 12.01         0.360"},
{"C5","C5 12.01         0.360"},{"C4","C4 12.01         0.360"},{"C0","C0 40.08              "},{"H ","H  1.008         0.161"},{"HC","HC 1.008         0.135"},
{"H1","H1 1.008         0.135"},{"H2","H2 1.008         0.135"},{"H3","H3 1.008         0.135"},{"HA","HA 1.008         0.167"},{"H4","H4 1.008         0.167"},
{"H5","H5 1.008         0.167"},{"HO","HO 1.008         0.135"},{"HS","HS 1.008         0.135"},{"HW","HW 1.008         0.000"},{"HP","HP 1.008         0.135"},
{"HZ","HZ 1.008         0.161"},{"F ","F  19.00         0.320"},{"Cl","Cl 35.45         1.910"},{"Br","Br 79.90         2.880"},{"I ","I  126.9         4.690"},
{"MG","MG 24.305        0.120"},{"N ","N  14.01         0.530"},{"NA","NA 14.01         0.530"},{"NB","NB 14.01         0.530"},{"NC","NC 14.01         0.530"},
{"N2","N2 14.01         0.530"},{"N3","N3 14.01         0.530"},{"NT","NT 14.01         0.530"},{"N*","N* 14.01         0.530"},{"NY","NY 14.01         0.530"},
{"O ","O  16.00         0.434"},{"O2","O2 16.00         0.434"},{"OW","OW 16.00         0.000"},{"OH","OH 16.00         0.465"},{"OS","OS 16.00         0.465"},
{"OP","OP 16.00         0.465"},{"P ","P  30.97         1.538"},{"S ","S  32.06         2.900"},{"SH","SH 32.06         2.900"},{"SO","SO 32.06"},
{"SS","SS 32.06"},{"CU","CU 63.55              "},{"FE","FE 55.00              "},{"Se","Se 78.96"},{"SG","SG 32.06         2.900"},{"N1","N1 14.01         0.530"},};

std::vector <std::string> mass_library_keys={
{"C "},{"CA"},{"CB"},{"CC"},{"CD"},{"CI"},{"CK"},{"CP"},{"CM"},{"CS"},
{"CN"},{"CQ"},{"CR"},{"CT"},{"CV"},{"CW"},{"C*"},{"CX"},{"CY"},{"CZ"},
{"C5"},{"C4"},{"C0"},{"H "},{"HC"},{"H1"},{"H2"},{"H3"},{"HA"},{"H4"},
{"H5"},{"HO"},{"HS"},{"HW"},{"HP"},{"HZ"},{"F "},{"Cl"},{"Br"},{"I "},
{"MG"},{"N "},{"NA"},{"NB"},{"NC"},{"N2"},{"N3"},{"NT"},{"N*"},{"NY"},
{"O "},{"O2"},{"OW"},{"OH"},{"OS"},{"OP"},{"P "},{"S "},{"SH"},{"SO"},
{"SS"},{"CU"},{"FE"},{"Se"},{"SG"},{"N1"}};


// Bond library (map of bond-pair strings to full frcmod string)
std::map<std::string, std::string> bond_defaults = {
{"H-O_sb"     ,"  553.00   0.970"}, {"C_sp2-H"   ,"  367.00   1.120"}, {"C_sp3-H"    ,"  340.00   1.120"}, {"C_sp-H"     ,"  400.00   1.120"}, {"H-N_sp2"    ,"  434.00   1.010"}, 
{"H-N_sp3"    ,"  434.00   1.010"}, {"H-S"       ,"  274.00   1.341"}, {"H-H"        ,"  553.00   0.737"}, {"C_sp2-C_sp2","  410.00   1.340"}, {"C_sp2-C_sp3","  317.00   1.500"}, 
{"C_sp2-N_sp2","  454.00   1.279"}, {"C_sp2-O_db","  613.00   1.230"}, {"C_sp2-O_sb" ,"  450.00   1.360"}, {"C_sp2-F"    ,"  386.00   1.380"}, {"C_sp2-Cl"   ,"  193.00   1.766"}, 
{"C_sp2-I"    ,"  171.00   2.166"}, {"Br-C_sp2"  ,"  172.00   1.944"}, {"C_sp3-C_sp3","  310.00   1.540"}, {"C_sp3-N_sp2","  337.00   1.465"}, {"C_sp3-O_sb" ,"  320.00   1.470"},
{"C_sp3-N_sp3","  367.00   1.469"}, {"C_sp3-S"   ,"  232.00   1.810"}, {"C_sp3-F"    ,"  367.00   1.380"}, {"C_sp3-Cl"   ,"  232.00   1.766"}, {"C_sp3-I"    ,"  148.00   2.166"}, 
{"Br-C_sp3"   ,"  159.00   1.944"}, {"C_sp-C_sp3","  400.00   1.460"}, {"C_sp-N_sp"  ,"  600.00   1.154"}, {"C_sp-C_sp"  ,"  600.00   1.206"}, {"O_db-P"     ,"  525.00   1.530"}, 
{"O_sb-P"     ,"  230.00   1.530"}, {"N_sp2-P"   ,"  250.00   1.840"}, {"S-S"        ,"  166.00   2.038"}, {"EP-O_db"    ,"  600.00   0.200"}, {"EP-O_sb"    ,"  600.00   0.200"},
{"EP-N_sp3"   ,"  600.00   0.200"}, {"EP-N_sp2"  ,"  600.00   0.200"}, {"EP-S"       ,"  600.00   0.700"}};

std::vector<std::string> bond_defaults_keys = {
"H-O_sb"      ,"C_sp2-H"    ,"C_sp3-H"    ,"C_sp-H"      ,"H-N_sp2",
"H-N_sp3"     ,"H-S"        ,"H-H"        ,"C_sp2-C_sp2" ,"C_sp2-C_sp3",
"C_sp2-N_sp2" ,"C_sp2-O_db" ,"C_sp2-O_sb" ,"C_sp2-F"     ,"C_sp2-Cl",
"C_sp2-I"     ,"Br-C_sp2"   ,"C_sp3-C_sp3","C_sp3-N_sp2" ,"C_sp3-O_sb",
"C_sp3-N_sp3" ,"C_sp3-S"    ,"C_sp3-F"    ,"C_sp3-Cl"    ,"C_sp3-I",
"Br-C_sp3"    ,"C_sp-C_sp3" ,"C_sp-N_sp"  ,"C_sp-C_sp"   ,"O_db-P",
"O_sb-P"      ,"N_sp2-P"    ,"S-S"        ,"EP-O_db"     ,"EP-O_sb",
"EP-N_sp3"    ,"EP-N_sp2"   ,"EP-S"};

std::map <std::string,std::string> bond_library={
{"C -C*","C -C*  447.00   1.419"},{"C -CD","C -CD  469.00   1.400"},{"C -CK","C -CK  450.00   1.400"},{"C -CQ","C -CQ  450.00   1.400"},{"C -CZ","C -CZ  450.00   1.400"},
{"C -H1","C -H1  367.00   1.080"},{"C -H2","C -H2  367.00   1.080"},{"C -HA","C -HA  367.00   1.080"},{"C -HC","C -HC  367.00   1.080"},{"C -N2","C -N2  450.00   1.400"},
{"C -NB","C -NB  418.00   1.388"},{"C -NT","C -NT  450.00   1.400"},{"CA-N ","CA-N   490.00   1.335"},{"CA-NB","CA-NB  418.00   1.388"},{"CA-NB","CA-NB  490.0    1.335"},
{"CB-N ","CB-N   418.00   1.388"},{"CB-N ","CB-N   490.0    1.335"},{"CC-N ","CC-N   490.00   1.335"},{"CD-N ","CD-N   450.00   1.400"},{"CD-NB","CD-NB  450.00   1.400"},
{"CK-N ","CK-N   490.0    1.335"},{"CM-N ","CM-N   440.00   1.371"},{"CM-NB","CM-NB  448.00   1.365"},{"CQ-N ","CQ-N   369.10   1.369"},{"CQ-NB","CQ-NB  450.00   1.400"},
{"CR-N ","CR-N   418.00   1.388"},{"CT-NB","CT-NB  436.00   1.374"},{"CZ-N ","CZ-N   424.00   1.383"},{"CZ-NB","CZ-NB  450.00   1.400"},{"H -NB","H -NB  434.00   1.010"},
{"N -OS","N -OS  450.00   1.400"},{"N*-NB","N*-NB  450.00   1.400"},{"N2-NB","N2-NB  450.00   1.400"},{"NB-NT","NB-NT  450.00   1.400"}};

std::vector <std::string> bond_library_keys={
{"C -C*"}, {"C -CD"}, {"C -CK"}, {"C -CQ"}, {"C -CZ"}, {"C -H1"}, {"C -H2"}, {"C -HA"}, {"C -HC"}, {"C -N2"},
{"C -NB"}, {"C -NT"}, {"CA-N "}, {"CA-NB"}, {"CA-NB"}, {"CB-N "}, {"CB-N "}, {"CC-N "}, {"CD-N "}, {"CD-NB"},
{"CK-N "}, {"CM-N "}, {"CM-NB"}, {"CQ-N "}, {"CQ-NB"}, {"CR-N "}, {"CT-NB"}, {"CZ-N "}, {"CZ-NB"}, {"H -NB"},
{"N -OS"}, {"N*-NB"}, {"N2-NB"}, {"NB-NT"}};


// Angle library (map of bond-pair strings to full frcmod string)
//      FC          DEG
// "   80.000     120.000"
std::map<std::string, std::string> angle_defaults = {
//sp carbon centers
{"C_sp3-C_sp-N_sp"  ,"   80.000     180.000"}, {"C_sp-C_sp-C_sp3"  ,"   80.000     180.000"}, {"C_sp-C_sp-H"      ,"   50.000     180.000"},
//sp2 carbon centers
{"C_sp2-C_sp2-O_db" ,"   80.000     124.900"}, {"C_sp2-C_sp2-O_sb" ,"   75.000     121.700"}, {"C_sp2-C_sp2-C_sp2","   63.000     118.900"}, {"C_sp2-C_sp2-N_sp2","   70.000     119.100"},
{"C_sp3-C_sp2-O_db" ,"   75.000     118.700"}, {"C_sp3-C_sp2-N_sp2","   70.000     118.300"}, {"C_sp3-C_sp2-C_sp3","   63.000     117.000"}, {"C_sp3-C_sp2-O_sb" ,"   68.000     110.800"},
{"N_sp2-C_sp2-N_sp2","   70.000     119.200"}, {"N_sp2-C_sp2-O_db" ,"   80.000     121.700"}, {"O_db-C_sp2-O_db"  ,"   80.000     126.000"}, {"O_db-C_sp2-O_sb"  ,"   80.000     122.500"},
{"C_sp2-C_sp2-H"    ,"   50.000     120.200"}, {"C_sp3-C_sp2-H"    ,"   50.000     116.100"}, {"H-C_sp2-O_db"     ,"   50.000     119.500"}, {"H-C_sp2-O_sb"     ,"   50.000     111.300"},
{"H-C_sp2-N_sp2"    ,"   50.000     120.700"}, {"C_sp2-C_sp2-C_sp3","   70.000     120.800"}, {"C_sp2-C_sp2-F"    ,"   70.000     121.000"}, {"C_sp2-C_sp2-Cl"   ,"   70.000     118.800"},
{"Br-C_sp2-C_sp2"   ,"   70.000     118.800"}, {"C_sp2-C_sp2-I"    ,"   70.000     118.800"}, {"H-C_sp2-H"        ,"   35.000     117.100"},
//sp3 carbon centers
{"H-C_sp3-H"        ,"   35.000     109.500"}, {"H-C_sp3-N_sp2"    ,"   50.000     109.500"}, {"H-C_sp3-O_sb"     ,"   50.000     109.500"}, {"C_sp2-C_sp3-H"    ,"   50.000     109.500"},
{"C_sp-C_sp3-H"     ,"   50.000     110.000"}, {"H-C_sp3-S"        ,"   50.000     109.500"}, {"H-C_sp3-N_sp3"    ,"   50.000     109.500"}, {"C_sp2-C_sp3-N_sp2","   63.000     110.100"},
{"C_sp2-C_sp3-N_sp3","   80.000     111.200"}, {"C_sp2-C_sp3-C_sp3","   61.100     112.500"}, {"C_sp2-C_sp3-O_sb" ,"   52.500     109.500"}, {"C_sp3-C_sp3-C_sp3","   40.000     109.500"},
{"C_sp3-C_sp3-H"    ,"   50.000     109.500"}, {"C_sp3-C_sp3-N_sp2","   72.500     110.000"}, {"C_sp3-C_sp3-O_sb" ,"   50.000     109.500"}, {"C_sp3-C_sp3-S"    ,"   50.000     111.700"},
{"C_sp3-C_sp3-N_sp3","   80.000     111.200"}, {"C_sp-C_sp3-C_sp3" ,"   63.000     110.000"}, {"O_sb-C_sp3-O_sb"  ,"  160.000     101.000"}, {"C_sp-C_sp3-O_sb"  ,"   50.000     110.000"},
{"N_sp2-C_sp3-O_sb" ,"   50.000     109.500"}, {"F-C_sp3-F"        ,"   77.000     109.100"}, {"F-C_sp3-H"        ,"   50.000     109.500"}, {"C_sp3-C_sp3-F"    ,"   50.000     109.000"},
{"C_sp3-C_sp3-Cl"   ,"   50.000     108.500"}, {"Cl-C_sp3-H"       ,"   50.000     108.500"}, {"Br-C_sp3-C_sp3"   ,"   50.000     108.000"}, {"Br-C_sp3-H"       ,"   50.000     106.500"},
{"C_sp3-C_sp3-I"    ,"   50.000     106.000"},
//sp2 nitrogen centers
{"C_sp2-N_sp2-P"    ,"   76.700     125.100"}, {"C_sp2-N_sp2-C_sp3","   64.500     123.700"}, {"C_sp2-N_sp2-H"    ,"   50.000     121.600"}, {"C_sp3-N_sp2-H"    ,"   50.000     118.200"},
{"C_sp3-N_sp2-C_sp3","   50.000     118.000"}, {"H-N_sp2-H"        ,"   35.000     120.000"}, {"C_sp2-N_sp2-C_sp2","   70.000     114.600"}, {"C_sp2-N_sp2-EP"   ,"  150.000     123.600"},
//sp3 nitrogen centers
{"C_sp3-N_sp3-H"    ,"   50.000     109.500"}, {"C_sp3-N_sp3-C_sp3","   50.000     109.500"}, {"H-N_sp3-H"        ,"   35.000     109.500"},
{"C_sp3-N_sp3-EP"   ,"  150.000     109.500"}, {"EP-N_sp3-H"       ,"  150.000     109.500"},
//single bond oxygen centers
{"H-O_sb-H"         ,"  100.000     104.500"}, {"C_sp2-O_sb-H"     ,"   51.700     111.500"}, {"C_sp3-O_sb-H"     ,"   55.000     108.500"}, {"H-O_sb-P"         ,"   45.000     108.500"},
{"C_sp2-O_sb-C_sp3" ,"   60.000     117.000"}, {"C_sp3-O_sb-C_sp3" ,"   60.000     109.500"}, {"C_sp3-O_sb-P"     ,"  100.000     120.500"}, {"C_sp2-O_sb-P"     ,"  100.000     120.500"},
{"P-O_sb-P"         ,"  100.000     120.500"}, {"C_sp2-O_sb-EP"    ,"  150.000     112.100"}, {"C_sp3-O_sb-EP"    ,"  150.000     109.500"}, {"EP-O_sb-H"        ,"  150.000     109.500"}, 
{"EP-O_sb-EP"       ,"  150.000     109.500"}, {"EP-O_sb-P"        ,"  150.000     109.500"},
//phosphorus centers
{"N_sp2-P-O_db"     ,"   42.900     102.400"}, {"O_db-P-O_sb"      ,"   81.700     108.200"}, {"O_db-P-O_db"      ,"  140.000     119.900"}, {"O_sb-P-O_sb"      ,"   45.000     102.600"},
//sulfur centers
{"C_sp3-S-C_sp3"    ,"   62.000      98.900"}, {"C_sp3-S-S"        ,"   68.000     103.700"}, {"C_sp3-S-H"        ,"   43.000      96.000"}, {"H-S-H"            ,"   35.000      92.100"}, 
{"C_sp3-S-EP"       ,"  150.000      90.000"}, {"EP-S-EP"          ,"  150.000     180.000"}, {"EP-S-H"           ,"  150.000      90.000"}, {"EP-S-S"           ,"  150.000      96.700"},
//double bond oxygen centers
{"C_sp2-O_db-EP"    ,"  150.000     120.000"}, {"EP-O_db-EP"       ,"  150.000     120.000"},
//defaults (unmatched groups)
{"O_sb"             ,"  104.000     112.700"}, {"H"                ,"    0.000     127.700"}, {"C_sp2"            ,"   62.300     119.300"}, {"N_sp2"            ,"   74.700     119.900"},
{"P"                ,"   82.200     109.000"}, {"C_sp3"            ,"   54.800     109.500"}, {"C_sp"             ,"   70.000     180.000"}, {"N_sp3"            ,"   80.800     109.500"},
{"S"                ,"  110.800     111.700"}, {"O_db"             ,"  150.000     120.000"}, {"H-H-O_sb"         ,"    0.000     127.700"} };

std::vector<std::string> angle_defaults_keys = {
"C_sp3-C_sp-N_sp"  , "C_sp-C_sp-C_sp3"  , "C_sp-C_sp-H"      ,  
"C_sp2-C_sp2-O_db" , "C_sp2-C_sp2-O_sb" , "C_sp2-C_sp2-C_sp2", "C_sp2-C_sp2-N_sp2", 
"C_sp3-C_sp2-O_db" , "C_sp3-C_sp2-N_sp2", "C_sp3-C_sp2-C_sp3", "C_sp3-C_sp2-O_sb" , 
"N_sp2-C_sp2-N_sp2", "N_sp2-C_sp2-O_db" , "O_db-C_sp2-O_db"  , "O_db-C_sp2-O_sb"  , 
"C_sp2-C_sp2-H"    , "C_sp3-C_sp2-H"    , "H-C_sp2-O_db"     , "H-C_sp2-O_sb"     , 
"H-C_sp2-N_sp2"    , "C_sp2-C_sp2-C_sp3", "C_sp2-C_sp2-F"    , "C_sp2-C_sp2-Cl"   , 
"Br-C_sp2-C_sp2"   , "C_sp2-C_sp2-I"    , "H-C_sp2-H"        , 
"H-C_sp3-H"        , "H-C_sp3-N_sp2"    , "H-C_sp3-O_sb"     , "C_sp2-C_sp3-H"    , 
"C_sp-C_sp3-H"     , "H-C_sp3-S"        , "H-C_sp3-N_sp3"    , "C_sp2-C_sp3-N_sp2",
"C_sp2-C_sp3-N_sp3", "C_sp2-C_sp3-C_sp3", "C_sp2-C_sp3-O_sb" , "C_sp3-C_sp3-C_sp3", 
"C_sp3-C_sp3-H"    , "C_sp3-C_sp3-N_sp2", "C_sp3-C_sp3-O_sb" , "C_sp3-C_sp3-S"    , 
"C_sp3-C_sp3-N_sp3", "C_sp-C_sp3-C_sp3" , "O_sb-C_sp3-O_sb"  , "C_sp-C_sp3-O_sb"  ,
"N_sp2-C_sp3-O_sb" , "F-C_sp3-F"        , "F-C_sp3-H"        , "C_sp3-C_sp3-F"    ,
"C_sp3-C_sp3-Cl"   , "Cl-C_sp3-H"       , "Br-C_sp3-C_sp3"   , "Br-C_sp3-H"       ,
"C_sp3-C_sp3-I"    ,
"C_sp2-N_sp2-P"    , "C_sp2-N_sp2-C_sp3", "C_sp2-N_sp2-H"    , "C_sp3-N_sp2-H"    ,
"C_sp3-N_sp2-C_sp3", "H-N_sp2-H"        , "C_sp2-N_sp2-C_sp2", "C_sp2-N_sp2-EP"   ,
"C_sp3-N_sp3-H"    , "C_sp3-N_sp3-C_sp3", "H-N_sp3-H"        , 
"C_sp3-N_sp3-EP"   , "EP-N_sp3-H"       , 
"H-O_sb-H"         , "C_sp2-O_sb-H"     , "C_sp3-O_sb-H"     , "H-O_sb-P"         ,
"C_sp2-O_sb-C_sp3" , "C_sp3-O_sb-C_sp3" , "C_sp3-O_sb-P"     , "C_sp2-O_sb-P"     ,
"P-O_sb-P"         , "C_sp2-O_sb-EP"    , "C_sp3-O_sb-EP"    , "EP-O_sb-H"        ,
"EP-O_sb-EP"       , "EP-O_sb-P"        , 
"N_sp2-P-O_db"     , "O_db-P-O_sb"      , "O_db-P-O_db"      , "O_sb-P-O_sb"      ,
"C_sp3-S-C_sp3"    , "C_sp3-S-S"        , "C_sp3-S-H"        , "H-S-H"            ,
"C_sp3-S-EP"       , "EP-S-EP"          , "EP-S-H"           , "EP-S-S"           ,
"C_sp2-O_db-EP"    , "EP-O_db-EP"       , 
"O_sb"             , "H"                , "C_sp2"            , "N_sp2"            ,
"P"                , "C_sp3"            , "C_sp"             , "N_sp3"            ,
"S"                , "O_db"             , "H-H-O_sb" };



std::map <std::string,std::string> angle_library={
{"C -C -CA","C -C -CA   80.000     120.000"},{"C -C -CB","C -C -CB   63.000     120.000"},{"C -C -CD","C -C -CD   50.000     120.000"},{"C -C -CM","C -C -CM   50.000     120.000"},{"C -C -CT","C -C -CT   50.000     120.000"},
{"C -C -N*","C -C -N*   50.000     120.000"},{"C -C -NA","C -C -NA   70.000     111.300"},{"C -C -NB","C -C -NB   70.000     111.300"},{"C -C -OS","C -C -OS   70.000     120.000"},{"C*-C -H4","C*-C -H4   50.000     120.000"},
{"C*-C -OS","C*-C -OS   50.000     120.000"},{"CA-C -CM","CA-C -CM   50.000     120.000"},{"CA-C -CT","CA-C -CT   50.000     120.000"},{"CA-C -HA","CA-C -HA   50.000     120.000"},{"CA-C -N2","CA-C -N2   50.000     120.000"},
{"CA-C -NB","CA-C -NB   50.000     120.000"},{"C*-C -CB","C*-C -CB   50.000     120.000"},{"CA-C -CB","CA-C -CB   50.000     120.000"},{"CB-C -CT","CB-C -CT   50.000     120.000"},{"CB-C -N ","CB-C -N    63.000     120.000"},
{"CB-C -NB","CB-C -NB   70.000     111.300"},{"CD-C -CD","CD-C -CD   63.000     120.000"},{"CD-C -H2","CD-C -H2   50.000     120.000"},{"CD-C -HA","CD-C -HA   50.000     120.000"},{"CD-C -N ","CD-C -N    50.000     120.000"},
{"CD-C -OS","CD-C -OS   50.000     120.000"},{"CA-C -CK","CA-C -CK   50.000     120.000"},{"CD-C -CM","CD-C -CM   50.000     120.000"},{"CM-C -CM","CM-C -CM   50.000     120.000"},{"CM-C -CT","CM-C -CT   63.000     120.000"},
{"CM-C -CZ","CM-C -CZ   50.000     120.000"},{"CM-C -HA","CM-C -HA   50.000     120.000"},{"CM-C -HC","CM-C -HC   63.000     120.000"},{"CM-C -N ","CM-C -N    50.000     120.000"},{"CM-C -N2","CM-C -N2   50.000     120.000"},
{"CM-C -NB","CM-C -NB   70.000     114.100"},{"CN-C -N ","CN-C -N    50.000     120.000"},{"CA-C -CQ","CA-C -CQ   50.000     120.000"},{"CT-C -H1","CT-C -H1   63.000     120.000"},{"CT-C -HC","CT-C -HC   63.000     120.000"},
{"CT-C -NB","CT-C -NB   70.000     116.600"},{"CT-C -CZ","CT-C -CZ   50.000     120.000"},{"H1-C -H1","H1-C -H1   63.000     120.000"},{"HC-C -HC","HC-C -HC   63.000     120.000"},{"C -C -N ","C -C -N    50.000     120.000"},
{"CA-C -N ","CA-C -N    50.000     120.000"},{"H1-C -N ","H1-C -N    63.000     120.000"},{"H2-C -N ","H2-C -N    50.000     120.000"},{"HA-C -N ","HA-C -N    50.000     120.000"},{"N -C -N*","N -C -N*   70.000     118.600"},
{"N -C -NA","N -C -NA   70.000     118.600"},{"N -C -NB","N -C -NB   70.000     118.600"},{"N -C -NC","N -C -NC   70.000     120.000"},{"N -C -OS","N -C -OS   80.000     122.900"},{"CM-C -N*","CM-C -N*   50.000     120.000"},
{"N*-C -NB","N*-C -NB   50.000     120.000"},{"C*-C -NB","C*-C -NB   50.000     120.000"},{"CD-C -NB","CD-C -NB   50.000     120.000"},{"H1-C -NB","H1-C -NB   50.000     120.000"},{"CD-C -NC","CD-C -NC   50.000     120.000"},
{"CM-C -NC","CM-C -NC   50.000     120.000"},{"CT-C -NT","CT-C -NT   50.000     120.000"},{"H1-C -NT","H1-C -NT   50.000     120.000"},{"NB-C -O ","NB-C -O    80.000     120.600"},{"CM-C -OS","CM-C -OS   50.000     120.000"},
{"H2-C -OS","H2-C -OS   50.000     120.000"},{"H4-C -OS","H4-C -OS   50.000     120.000"},{"H1-C -OH","H1-C -OH   50.000     120.000"},{"CT-C -H2","CT-C -H2   63.000     120.000"},{"C -C*-C*","C -C*-C*   70.000     108.700"},
{"C -C*-CC","C -C*-CC   70.000     108.700"},{"C -C*-CM","C -C*-CM   70.000     108.700"},{"C -C*-CQ","C -C*-CQ   70.000     108.700"},{"C -C*-CT","C -C*-CT   70.000     108.700"},{"C -C*-HA","C -C*-HA   70.000     108.700"},
{"C -CA-C ","C -CA-C    80.000     120.000"},{"C -CA-CD","C -CA-CD   50.000     120.000"},{"C -CA-CM","C -CA-CM   80.000     120.000"},{"C -CA-CT","C -CA-CT   50.000     120.000"},{"CA-CA-N ","CA-CA-N    70.000     120.000"},
{"CA-CA-NB","CA-CA-NB   70.000     120.000"},{"C -CA-CB","C -CA-CB   50.000     120.000"},{"CB-CA-N ","CB-CA-N    70.000     117.300"},{"C -CA-CK","C -CA-CK   50.000     120.000"},{"CK-CA-NB","CK-CA-NB   50.000     120.000"},
{"CM-CA-N ","CM-CA-N    70.000     120.100"},{"CM-CA-N ","CM-CA-N    80.000     120.000"},{"CM-CA-NB","CM-CA-NB   50.000     120.000"},{"C -CA-N ","C -CA-N    50.000     120.000"},{"CD-CA-N ","CD-CA-N    50.000     120.000"},
{"CT-CA-N ","CT-CA-N    50.000     120.000"},{"CZ-CA-N ","CZ-CA-N    50.000     120.000"},{"H4-CA-N ","H4-CA-N    50.000     120.000"},{"H5-CA-N ","H5-CA-N    50.000     120.000"},{"N -CA-N*","N -CA-N*   70.000     120.100"},
{"N -CA-N2","N -CA-N2   70.000     119.300"},{"N -CA-NB","N -CA-NB   80.000     120.000"},{"N -CA-NC","N -CA-NC   70.000     123.300"},{"N -CA-OS","N -CA-OS   50.000     120.000"},{"C -CA-N*","C -CA-N*   50.000     120.000"},
{"C -CA-NA","C -CA-NA   50.000     120.000"},{"NA-CA-NB","NA-CA-NB   70.000     116.000"},{"C*-CA-NB","C*-CA-NB   50.000     120.000"},{"CD-CA-NB","CD-CA-NB   50.000     120.000"},{"H4-CA-NB","H4-CA-NB   50.000     120.000"},
{"N2-CA-NB","N2-CA-NB   70.000     116.000"},{"NB-CA-NC","NB-CA-NC   80.000     120.000"},{"C -CA-NC","C -CA-NC   50.000     120.000"},{"C -CA-OS","C -CA-OS   50.000     120.000"},{"N -CA-NA","N -CA-NA   80.000     120.000"},
{"C -CA-N2","C -CA-N2   70.000     130.000"},{"C -CB-C ","C -CB-C    50.000     120.000"},{"C -CB-CA","C -CB-CA   50.000     120.000"},{"C -CB-N ","C -CB-N    70.000     132.400"},{"CB-CB-N ","CB-CB-N    70.000     110.400"},
{"CM-CB-NB","CM-CB-NB   50.000     120.000"},{"N -CB-NB","N -CB-NB   70.000     120.000"},{"C -CB-OS","C -CB-OS   50.000     120.000"},{"N -CB-N*","N -CB-N*   70.000     120.000"},{"N -CB-NC","N -CB-NC   70.000     120.000"},
{"C*-CC-N ","C*-CC-N    70.000     108.700"},{"CA-CC-N ","CA-CC-N    70.000     108.700"},{"C -CD-CD","C -CD-CD   63.000     120.000"},{"C -CD-CM","C -CD-CM   50.000     120.000"},{"C -CD-HA","C -CD-HA   50.000     120.000"},
{"CD-CD-N ","CD-CD-N    50.000     120.000"},{"CD-CD-NB","CD-CD-NB   50.000     120.000"},{"CM-CD-NB","CM-CD-NB   50.000     120.000"},{"CA-CD-N ","CA-CD-N    50.000     120.000"},{"C -CK-N*","C -CK-N*   50.000     120.000"},
{"C -CK-NB","C -CK-NB   50.000     120.000"},{"NB-CK-NB","NB-CK-NB   70.000     113.900"},{"H5-CK-N ","H5-CK-N    50.000     120.000"},{"N -CK-NB","N -CK-NB   50.000     120.000"},{"C -CM-C ","C -CM-C    50.000     120.000"},
{"C -CM-CA","C -CM-CA   50.000     120.000"},{"C -CM-CD","C -CM-CD   50.000     120.000"},{"C -CM-CQ","C -CM-CQ   50.000     120.000"},{"C -CM-N ","C -CM-N    63.000     117.000"},{"C -CM-NB","C -CM-NB   63.000     117.000"},
{"CA-CM-N ","CA-CM-N    63.000     117.000"},{"CA-CM-NB","CA-CM-NB   63.000     117.000"},{"CB-CM-N ","CB-CM-N    50.000     120.000"},{"CC-CM-N ","CC-CM-N    63.000     117.000"},{"CD-CM-N ","CD-CM-N    50.000     120.000"},
{"CD-CM-NB","CD-CM-NB   50.000     120.000"},{"C -CM-CK","C -CM-CK   50.000     120.000"},{"CM-CM-N ","CM-CM-N    63.000     117.000"},{"CM-CM-NB","CM-CM-NB   63.000     117.000"},{"CQ-CM-NB","CQ-CM-NB   50.000     120.000"},
{"CT-CM-N ","CT-CM-N    50.000     120.000"},{"CT-CM-NB","CT-CM-NB   63.000     117.000"},{"CZ-CM-N ","CZ-CM-N    50.000     120.000"},{"H4-CM-N ","H4-CM-N    70.000     119.300"},{"HA-CM-N ","HA-CM-N    50.000     120.000"},
{"N -CM-N ","N -CM-N    50.000     120.000"},{"N -CM-N2","N -CM-N2   50.000     120.000"},{"N -CM-NA","N -CM-NA   63.000     117.000"},{"N -CM-NB","N -CM-NB   50.000     120.000"},{"N -CM-NT","N -CM-NT   63.000     117.000"},
{"N -CM-OS","N -CM-OS   50.000     120.000"},{"C -CM-N*","C -CM-N*   50.000     120.000"},{"N2-CM-NB","N2-CM-NB   50.000     120.000"},{"C -CM-NA","C -CM-NA   50.000     120.000"},{"C -CM-NC","C -CM-NC   50.000     120.000"},
{"N -CM-OH","N -CM-OH   50.000     120.000"},{"C -CM-OS","C -CM-OS   50.000     120.000"},{"C*-CQ-N ","C*-CQ-N    70.000     108.700"},{"CM-CQ-N ","CM-CQ-N    70.000     108.700"},{"H5-CQ-N ","H5-CQ-N    50.000     120.000"},
{"C -CQ-N ","C -CQ-N    50.000     120.000"},{"CD-CQ-N ","CD-CQ-N    50.000     120.000"},{"CT-CQ-N ","CT-CQ-N    50.000     120.000"},{"N -CQ-N2","N -CQ-N2   50.000     120.000"},{"N -CQ-OS","N -CQ-OS   50.000     120.000"},
{"CA-CQ-NB","CA-CQ-NB   50.000     120.000"},{"CD-CQ-NB","CD-CQ-NB   50.000     120.000"},{"C -CQ-NC","C -CQ-NC   50.000     120.000"},{"N -CQ-NC","N -CQ-NC   70.000     129.100"},{"NB-CQ-NC","NB-CQ-NC   50.000     120.000"},
{"H5-CR-N ","H5-CR-N    50.000     120.000"},{"N -CR-NB","N -CR-NB   70.000     120.000"},{"C -CT-C ","C -CT-C    50.000     109.500"},{"C -CT-CA","C -CT-CA   50.000     109.500"},{"C -CT-CD","C -CT-CD   50.000     109.500"},
{"C -CT-CM","C -CT-CM   50.000     109.500"},{"C -CT-CZ","C -CT-CZ   50.000     109.500"},{"C -CT-NB","C -CT-NB   50.000     109.500"},{"C -CT-NT","C -CT-NT   50.000     109.500"},{"C -CT-O ","C -CT-O    50.000     109.500"},
{"C -CT-OH","C -CT-OH   50.000     109.500"},{"CA-CT-NB","CA-CT-NB   50.000     109.500"},{"C -CT-CB","C -CT-CB   50.000     109.500"},{"CB-CT-N ","CB-CT-N    50.000     109.500"},{"C -CT-CK","C -CT-CK   50.000     109.500"},
{"CK-CT-N ","CK-CT-N    70.000     109.500"},{"CK-CT-NB","CK-CT-NB   50.000     109.500"},{"CM-CT-NB","CM-CT-NB   50.000     109.500"},{"CQ-CT-NB","CQ-CT-NB   50.000     109.500"},{"CT-CT-NB","CT-CT-NB   50.000     109.500"},
{"CA-CT-N ","CA-CT-N    50.000     109.500"},{"CM-CT-N ","CM-CT-N    50.000     109.500"},{"H2-CT-N ","H2-CT-N    70.000     109.500"},{"HP-CT-N ","HP-CT-N    50.000     109.500"},{"N -CT-N ","N -CT-N    63.000     110.100"},
{"N -CT-NB","N -CT-NB   50.000     109.500"},{"N -CT-NT","N -CT-NT   50.000     109.500"},{"N -CT-OS","N -CT-OS   70.000     109.500"},{"C -CT-N*","C -CT-N*   50.000     109.500"},{"H1-CT-NB","H1-CT-NB   50.000     109.500"},
{"NB-CT-P ","NB-CT-P    50.000     109.500"},{"C -CT-H2","C -CT-H2   50.000     109.500"},{"C -CZ-NB","C -CZ-NB   50.000     120.000"},{"N -CZ-NB","N -CZ-NB   50.000     120.000"},{"C -N -C ","C -N -C    50.000     120.000"},
{"C -N -CA","C -N -CA   50.000     121.900"},{"C -N -CB","C -N -CB   50.000     121.900"},{"C -N -CM","C -N -CM   50.000     120.000"},{"CA-N -CA","CA-N -CA   50.000     121.900"},{"CA-N -CM","CA-N -CM   70.000     120.500"},
{"CA-N -CQ","CA-N -CQ   70.000     118.600"},{"CA-N -CT","CA-N -CT   50.000     121.900"},{"CA-N -H ","CA-N -H    50.000     119.100"},{"CA-N -OS","CA-N -OS   50.000     120.000"},{"CB-N -CR","CB-N -CR   50.000     118.000"},
{"CB-N -CT","CB-N -CT   50.000     121.900"},{"CB-N -H ","CB-N -H    50.000     120.000"},{"CC-N -CM","CC-N -CM   50.000     120.000"},{"CD-N -CT","CD-N -CT   50.000     120.000"},{"CM-N -CM","CM-N -CM   50.000     120.000"},
{"CM-N -CT","CM-N -CT   70.000     118.600"},{"CM-N -H ","CM-N -H    50.000     120.000"},{"CM-N -OS","CM-N -OS   50.000     120.000"},{"CM-N -CQ","CM-N -CQ   50.000     120.000"},{"CR-N -H ","CR-N -H    50.000     120.000"},
{"CT-N -OS","CT-N -OS   50.000     120.000"},{"CT-N -CZ","CT-N -CZ   50.000     120.000"},{"CZ-N -H ","CZ-N -H    50.000     120.000"},{"CB-N -CK","CB-N -CK   50.000     120.000"},{"CK-N -CT","CK-N -CT   50.000     120.000"},
{"CA-N -CB","CA-N -CB   50.000     121.900"},{"C -N -CK","C -N -CK   50.000     120.000"},{"C -N*-CB","C -N*-CB   50.000     120.000"},{"C -N*-CR","C -N*-CR   50.000     120.000"},{"CA-N*-NB","CA-N*-NB   50.000     120.000"},
{"CM-N*-NB","CM-N*-NB   50.000     120.000"},{"C -N2-CT","C -N2-CT   50.000     120.000"},{"CA-N2-NB","CA-N2-NB   50.000     120.000"},{"CM-N2-NB","CM-N2-NB   50.000     120.000"},{"CQ-N2-NB","CQ-N2-NB   50.000     120.000"},
{"CR-N2-NB","CR-N2-NB   50.000     120.000"},{"H -N2-NB","H -N2-NB   50.000     120.000"},{"C -NB-C ","C -NB-C    70.000     120.500"},{"C -NB-CR","C -NB-CR   70.000     117.000"},{"C -NB-CT","C -NB-CT   50.000     120.000"},
{"C -NB-H ","C -NB-H    70.000     120.500"},{"C -NB-N2","C -NB-N2   50.000     120.000"},{"C -NB-NT","C -NB-NT   50.000     120.000"},{"C -NB-CA","C -NB-CA   50.000     120.000"},{"CA-NB-CA","CA-NB-CA   50.000     120.000"},
{"CA-NB-CB","CA-NB-CB   70.000     117.000"},{"CA-NB-CM","CA-NB-CM   70.000     103.800"},{"CA-NB-CT","CA-NB-CT   50.000     120.000"},{"CA-NB-H ","CA-NB-H    70.000     120.500"},{"CB-NB-CT","CB-NB-CT   70.000     103.800"},
{"CD-NB-CT","CD-NB-CT   50.000     120.000"},{"CK-NB-CT","CK-NB-CT   70.000     103.800"},{"CB-NB-CM","CB-NB-CM   50.000     120.000"},{"CD-NB-CM","CD-NB-CM   50.000     120.000"},{"CM-NB-CM","CM-NB-CM   70.000     120.500"},
{"CM-NB-CT","CM-NB-CT   50.000     120.000"},{"CM-NB-H ","CM-NB-H    70.000     120.500"},{"C -NB-CQ","C -NB-CQ   50.000     120.000"},{"CT-NB-CT","CT-NB-CT   50.000     120.000"},{"C -NB-CZ","C -NB-CZ   50.000     120.000"},
{"CZ-NB-H ","CZ-NB-H    50.000     120.000"},{"H -NB-H ","H -NB-H    70.000     119.000"},{"H -NB-N2","H -NB-N2   50.000     120.000"},{"H -NB-NT","H -NB-NT   50.000     120.000"},{"C -NB-N*","C -NB-N*   50.000     120.000"},
{"CZ-NB-N*","CZ-NB-N*   50.000     120.000"},{"C -NC-CM","C -NC-CM   50.000     120.000"},{"C -NT-CT","C -NT-CT   50.000     109.500"},{"H -NT-NB","H -NT-NB   50.000     109.500"},{"C -OS-CA","C -OS-CA   63.000     117.000"},
{"C -OS-CM","C -OS-CM   63.000     117.000"},{"C -OS-N ","C -OS-N    50.000     109.500"},{"N -OS-SO","N -OS-SO   50.000     109.500"},};

std::vector <std::string> angle_library_keys={
{"C -C -CA"},{"C -C -CB"},{"C -C -CD"},{"C -C -CM"},{"C -C -CT"},{"C -C -N*"},{"C -C -NA"},{"C -C -NB"},{"C -C -OS"},{"C*-C -H4"},
{"C*-C -OS"},{"CA-C -CM"},{"CA-C -CT"},{"CA-C -HA"},{"CA-C -N2"},{"CA-C -NB"},{"C*-C -CB"},{"CA-C -CB"},{"CB-C -CT"},{"CB-C -N "},
{"CB-C -NB"},{"CD-C -CD"},{"CD-C -H2"},{"CD-C -HA"},{"CD-C -N "},{"CD-C -OS"},{"CA-C -CK"},{"CD-C -CM"},{"CM-C -CM"},{"CM-C -CT"},
{"CM-C -CZ"},{"CM-C -HA"},{"CM-C -HC"},{"CM-C -N "},{"CM-C -N2"},{"CM-C -NB"},{"CN-C -N "},{"CA-C -CQ"},{"CT-C -H1"},{"CT-C -HC"},
{"CT-C -NB"},{"CT-C -CZ"},{"H1-C -H1"},{"HC-C -HC"},{"C -C -N "},{"CA-C -N "},{"H1-C -N "},{"H2-C -N "},{"HA-C -N "},{"N -C -N*"},
{"N -C -NA"},{"N -C -NB"},{"N -C -NC"},{"N -C -OS"},{"CM-C -N*"},{"N*-C -NB"},{"C*-C -NB"},{"CD-C -NB"},{"H1-C -NB"},{"CD-C -NC"},
{"CM-C -NC"},{"CT-C -NT"},{"H1-C -NT"},{"NB-C -O "},{"CM-C -OS"},{"H2-C -OS"},{"H4-C -OS"},{"H1-C -OH"},{"CT-C -H2"},{"C -C*-C*"},
{"C -C*-CC"},{"C -C*-CM"},{"C -C*-CQ"},{"C -C*-CT"},{"C -C*-HA"},{"C -CA-C "},{"C -CA-CD"},{"C -CA-CM"},{"C -CA-CT"},{"CA-CA-N "},
{"CA-CA-NB"},{"C -CA-CB"},{"CB-CA-N "},{"C -CA-CK"},{"CK-CA-NB"},{"CM-CA-N "},{"CM-CA-N "},{"CM-CA-NB"},{"C -CA-N "},{"CD-CA-N "},
{"CT-CA-N "},{"CZ-CA-N "},{"H4-CA-N "},{"H5-CA-N "},{"N -CA-N*"},{"N -CA-N2"},{"N -CA-NB"},{"N -CA-NC"},{"N -CA-OS"},{"C -CA-N*"},
{"C -CA-NA"},{"NA-CA-NB"},{"C*-CA-NB"},{"CD-CA-NB"},{"H4-CA-NB"},{"N2-CA-NB"},{"NB-CA-NC"},{"C -CA-NC"},{"C -CA-OS"},{"N -CA-NA"},
{"C -CA-N2"},{"C -CB-C "},{"C -CB-CA"},{"C -CB-N "},{"CB-CB-N "},{"CM-CB-NB"},{"N -CB-NB"},{"C -CB-OS"},{"N -CB-N*"},{"N -CB-NC"},
{"C*-CC-N "},{"CA-CC-N "},{"C -CD-CD"},{"C -CD-CM"},{"C -CD-HA"},{"CD-CD-N "},{"CD-CD-NB"},{"CM-CD-NB"},{"CA-CD-N "},{"C -CK-N*"},
{"C -CK-NB"},{"NB-CK-NB"},{"H5-CK-N "},{"N -CK-NB"},{"C -CM-C "},{"C -CM-CA"},{"C -CM-CD"},{"C -CM-CQ"},{"C -CM-N "},{"C -CM-NB"},
{"CA-CM-N "},{"CA-CM-NB"},{"CB-CM-N "},{"CC-CM-N "},{"CD-CM-N "},{"CD-CM-NB"},{"C -CM-CK"},{"CM-CM-N "},{"CM-CM-NB"},{"CQ-CM-NB"},
{"CT-CM-N "},{"CT-CM-NB"},{"CZ-CM-N "},{"H4-CM-N "},{"HA-CM-N "},{"N -CM-N "},{"N -CM-N2"},{"N -CM-NA"},{"N -CM-NB"},{"N -CM-NT"},
{"N -CM-OS"},{"C -CM-N*"},{"N2-CM-NB"},{"C -CM-NA"},{"C -CM-NC"},{"N -CM-OH"},{"C -CM-OS"},{"C*-CQ-N "},{"CM-CQ-N "},{"H5-CQ-N "},
{"C -CQ-N "},{"CD-CQ-N "},{"CT-CQ-N "},{"N -CQ-N2"},{"N -CQ-OS"},{"CA-CQ-NB"},{"CD-CQ-NB"},{"C -CQ-NC"},{"N -CQ-NC"},{"NB-CQ-NC"},
{"H5-CR-N "},{"N -CR-NB"},{"C -CT-C "},{"C -CT-CA"},{"C -CT-CD"},{"C -CT-CM"},{"C -CT-CZ"},{"C -CT-NB"},{"C -CT-NT"},{"C -CT-O "},
{"C -CT-OH"},{"CA-CT-NB"},{"C -CT-CB"},{"CB-CT-N "},{"C -CT-CK"},{"CK-CT-N "},{"CK-CT-NB"},{"CM-CT-NB"},{"CQ-CT-NB"},{"CT-CT-NB"},
{"CA-CT-N "},{"CM-CT-N "},{"H2-CT-N "},{"HP-CT-N "},{"N -CT-N "},{"N -CT-NB"},{"N -CT-NT"},{"N -CT-OS"},{"C -CT-N*"},{"H1-CT-NB"},
{"NB-CT-P "},{"C -CT-H2"},{"C -CZ-NB"},{"N -CZ-NB"},{"C -N -C "},{"C -N -CA"},{"C -N -CB"},{"C -N -CM"},{"CA-N -CA"},{"CA-N -CM"},
{"CA-N -CQ"},{"CA-N -CT"},{"CA-N -H "},{"CA-N -OS"},{"CB-N -CR"},{"CB-N -CT"},{"CB-N -H "},{"CC-N -CM"},{"CD-N -CT"},{"CM-N -CM"},
{"CM-N -CT"},{"CM-N -H "},{"CM-N -OS"},{"CM-N -CQ"},{"CR-N -H "},{"CT-N -OS"},{"CT-N -CZ"},{"CZ-N -H "},{"CB-N -CK"},{"CK-N -CT"},
{"CA-N -CB"},{"C -N -CK"},{"C -N*-CB"},{"C -N*-CR"},{"CA-N*-NB"},{"CM-N*-NB"},{"C -N2-CT"},{"CA-N2-NB"},{"CM-N2-NB"},{"CQ-N2-NB"},
{"CR-N2-NB"},{"H -N2-NB"},{"C -NB-C "},{"C -NB-CR"},{"C -NB-CT"},{"C -NB-H "},{"C -NB-N2"},{"C -NB-NT"},{"C -NB-CA"},{"CA-NB-CA"},
{"CA-NB-CB"},{"CA-NB-CM"},{"CA-NB-CT"},{"CA-NB-H "},{"CB-NB-CT"},{"CD-NB-CT"},{"CK-NB-CT"},{"CB-NB-CM"},{"CD-NB-CM"},{"CM-NB-CM"},
{"CM-NB-CT"},{"CM-NB-H "},{"C -NB-CQ"},{"CT-NB-CT"},{"C -NB-CZ"},{"CZ-NB-H "},{"H -NB-H "},{"H -NB-N2"},{"H -NB-NT"},{"C -NB-N*"},
{"CZ-NB-N*"},{"C -NC-CM"},{"C -NT-CT"},{"H -NT-NB"},{"C -OS-CA"},{"C -OS-CM"},{"C -OS-N "},{"N -OS-SO"},};

// // Torsion library (map of bond-pair strings to full frcmod string)
// std::map <std::string, std::string> torsion_library;
// std::vector <std::string> torsion_library_keys;

// Dihedral library (map of bond-pair strings to full frcmod string)
std::map <std::string,std::string> DIHEDRAL_SET_DICT={
{"C_sp2-C_sp2-C_sp2-O_db","WW-XX-YY-ZZ   1    2.175       180.0            -2.\nWW-XX-YY-ZZ   1    0.30          0.0             3."},
{"C_sp3-C_sp2-C_sp2-C_sp3","WW-XX-YY-ZZ   1    6.65        180.0            -2.\nWW-XX-YY-ZZ   1    1.90        180.0             1."},
{"C_sp2-C_sp2-C_sp2-X ","WW-XX-YY-ZZ   1    2.175       180.0            -2.\nWW-XX-YY-ZZ   1    0.30          0.0             3."},
{"O_db-C_sp2-C_sp2-X ","WW-XX-YY-ZZ   1    2.175       180.0            -2.\nWW-XX-YY-ZZ   1    0.30          0.0             3."},
{"C_sp3-C_sp2-C_sp2-X ","WW-XX-YY-ZZ   1    6.65        180.0            -2.\nWW-XX-YY-ZZ   1    1.90        180.0             1."},
{"X -C_sp2-C_sp2-X ","WW-XX-YY-ZZ    4   16.10        180.0             2. "},
{"H -C_sp2-C_sp3-O_sb","WW-XX-YY-ZZ   1    0.00          0.0            -3.\nWW-XX-YY-ZZ   1    0.250         0.0             1."},
{"O_sb-C_sp2-C_sp3-C_sp3","WW-XX-YY-ZZ   1    1.178040    190.97653        -1. \nWW-XX-YY-ZZ   1    0.092102    295.63279        -2. \nWW-XX-YY-ZZ   1    0.962830    348.09535         3. "},
{"C_sp2-C_sp2-C_sp3-H ","WW-XX-YY-ZZ   1    0.38        180.0            -3. \nWW-XX-YY-ZZ   1    1.15          0.0             1. "},
{"N_sp2-C_sp2-C_sp3-N_sp2","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.55        180.0            -3. \nWW-XX-YY-ZZ   1    1.58        180.0            -2.\nWW-XX-YY-ZZ   1    0.45        180.0             1."},
{"N_sp2-C_sp2-C_sp3-C_sp3","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3. \nWW-XX-YY-ZZ   1    0.20          0.0            -2.\nWW-XX-YY-ZZ   1    0.20          0.0             1."},
{"O_db-C_sp2-C_sp3-H ","WW-XX-YY-ZZ   1    0.80          0.0            -1. \nWW-XX-YY-ZZ   1    0.00          0.0            -2. \nWW-XX-YY-ZZ   1    0.08        180.0             3. "},
{"X -C_sp2-C_sp3-O_sb","WW-XX-YY-ZZ   1    0.00          0.0            -3.\nWW-XX-YY-ZZ   1    0.250         0.0             1."},
{"X -C_sp2-C_sp3-C_sp3","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3. \nWW-XX-YY-ZZ   1    0.20          0.0            -2.\nWW-XX-YY-ZZ   1    0.20          0.0             1."},
{"X -C_sp2-C_sp3-H ","WW-XX-YY-ZZ   1    0.80          0.0            -1. \nWW-XX-YY-ZZ   1    0.00          0.0            -2. \nWW-XX-YY-ZZ   1    0.08        180.0             3. "},
{"X -C_sp2-C_sp3-N_sp2","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.55        180.0            -3. \nWW-XX-YY-ZZ   1    1.58        180.0            -2.\nWW-XX-YY-ZZ   1    0.45        180.0             1."},
{"H -C_sp2-C_sp3-X ","WW-XX-YY-ZZ   1    0.00          0.0            -3.\nWW-XX-YY-ZZ   1    0.250         0.0             1."},
{"O_sb-C_sp2-C_sp3-X ","WW-XX-YY-ZZ   1    1.178040    190.97653        -1. \nWW-XX-YY-ZZ   1    0.092102    295.63279        -2. \nWW-XX-YY-ZZ   1    0.962830    348.09535         3. "},
{"C_sp2-C_sp2-C_sp3-X ","WW-XX-YY-ZZ   1    0.38        180.0            -3. \nWW-XX-YY-ZZ   1    1.15          0.0             1. "},
{"N_sp2-C_sp2-C_sp3-X ","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.55        180.0            -3. \nWW-XX-YY-ZZ   1    1.58        180.0            -2.\nWW-XX-YY-ZZ   1    0.45        180.0             1."},
{"X -C_sp2-C_sp3-X ","WW-XX-YY-ZZ   6    0.00          0.0             2. "},
{"O_db-C_sp2-N_sp2-H ","WW-XX-YY-ZZ   1    2.50        180.0            -2. \nWW-XX-YY-ZZ   1    2.00          0.0             1. "},
{"O_db-C_sp2-N_sp2-EP","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"C_sp3-C_sp2-N_sp2-EP","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"X -C_sp2-N_sp2-H ","WW-XX-YY-ZZ   1    2.50        180.000          -2. \nWW-XX-YY-ZZ   1    2.00          0.0             1. "},
{"X -C_sp2-N_sp2-EP","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"O_db-C_sp2-N_sp2-X ","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"C_sp3-C_sp2-N_sp2-X ","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"X -C_sp2-N_sp2-X ","WW-XX-YY-ZZ   2    7.40        180.000           2. "},
{"X -C_sp2-O_db-X ","WW-XX-YY-ZZ   4   11.20        180.000           2. "},
{"O_db-C_sp2-O_sb-H ","WW-XX-YY-ZZ   1    2.30        180.000          -2.\nWW-XX-YY-ZZ   1    1.90          0.0             1."},
{"O_db-C_sp2-O_sb-C_sp3","WW-XX-YY-ZZ   1    2.70        180.000          -2.\nWW-XX-YY-ZZ   1    1.40        180.0             1."},
{"C_sp2-C_sp2-O_sb-EP","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"H -C_sp2-O_sb-EP","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"O_db-C_sp2-O_sb-X ","WW-XX-YY-ZZ   1    2.30        180.000          -2.\nWW-XX-YY-ZZ   1    1.90          0.0             1."},
{"C_sp2-C_sp2-O_sb-X ","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"H -C_sp2-O_sb-X ","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"X -C_sp2-O_sb-H ","WW-XX-YY-ZZ   1    2.30        180.000          -2.\nWW-XX-YY-ZZ   1    1.90          0.0             1."},
{"X -C_sp2-O_sb-C_sp3","WW-XX-YY-ZZ   1    2.70        180.000          -2.\nWW-XX-YY-ZZ   1    1.40        180.0             1."},
{"X -C_sp2-O_sb-EP","WW-XX-YY-ZZ   1    0.000       180.000           2.0"},
{"X -C_sp2-O_sb-X ","WW-XX-YY-ZZ   2    4.00        180.000           2."},
{"X -C_sp2-S-X ","WW-XX-YY-ZZ   2    4.00        180.000           2."},
{"X -C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."},
{"X -C_sp3-C_sp-X ","WW-XX-YY-ZZ   3    0.00          0.0             1. "},
{"O_sb-C_sp3-N_sp2-C_sp2","WW-XX-YY-ZZ   1    0.98023     110.0984         -1.\nWW-XX-YY-ZZ   1    1.38071      13.7765         -2.\nWW-XX-YY-ZZ   1    0.60481     176.3635         -3.\nWW-XX-YY-ZZ   1    0.30675      17.8068          4."},
{"C_sp2-C_sp3-N_sp2-C_sp2","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.42          0.0            -3.\nWW-XX-YY-ZZ   1    0.27          0.0            -2.\nWW-XX-YY-ZZ   1    0.00          0.0             1."},
{"C_sp3-C_sp3-N_sp2-C_sp2","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3.\nWW-XX-YY-ZZ   1    2.00          0.0            -2.\nWW-XX-YY-ZZ   1    2.00          0.0             1."},
{"X -C_sp3-N_sp2-C_sp2","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.42          0.0            -3.\nWW-XX-YY-ZZ   1    0.27          0.0            -2.\nWW-XX-YY-ZZ   1    0.00          0.0             1."},
{"O_sb-C_sp3-N_sp2-X ","WW-XX-YY-ZZ   1    0.98023     110.0984         -1.\nWW-XX-YY-ZZ   1    1.38071      13.7765         -2.\nWW-XX-YY-ZZ   1    0.60481     176.3635         -3.\nWW-XX-YY-ZZ   1    0.30675      17.8068          4."},
{"C_sp2-C_sp3-N_sp2-X ","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.42          0.0            -3.\nWW-XX-YY-ZZ   1    0.27          0.0            -2.\nWW-XX-YY-ZZ   1    0.00          0.0             1."},
{"C_sp3-C_sp3-N_sp2-X ","WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3.\nWW-XX-YY-ZZ   1    2.00          0.0            -2.\nWW-XX-YY-ZZ   1    2.00          0.0             1."},
{"X -C_sp3-N_sp2-X ","WW-XX-YY-ZZ   6    0.00          0.0             2. "},
{"C_sp3-C_sp3-N_sp3-C_sp3","WW-XX-YY-ZZ   1    0.30          0.0            -3.\nWW-XX-YY-ZZ   1    0.48        180.0             2."},
{"H -C_sp3-N_sp3-EP","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"C_sp3-C_sp3-N_sp3-EP","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"X -C_sp3-N_sp3-C_sp3","WW-XX-YY-ZZ   1    0.30          0.0            -3.\nWW-XX-YY-ZZ   1    0.48        180.0             2."},
{"X -C_sp3-N_sp3-EP","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"C_sp3-C_sp3-N_sp3-X ","WW-XX-YY-ZZ   1    0.30          0.0            -3.\nWW-XX-YY-ZZ   1    0.48        180.0             2."},
{"H -C_sp3-N_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"X -C_sp3-N_sp3-X ","WW-XX-YY-ZZ    6    1.80          0.0             3. "},
{"C_sp2-C_sp3-O_sb-C_sp3","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.100       180.0             2."},
{"C_sp3-C_sp3-O_sb-H ","WW-XX-YY-ZZ   1    0.16          0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1."},
{"C_sp3-C_sp3-O_sb-C_sp3","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.1         180.0             2."},
{"C_sp3-C_sp3-O_sb-C_sp2","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.80        180.0             1."},
{"C_sp3-C_sp3-O_sb-EP","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"O_sb-C_sp3-O_sb-C_sp3","WW-XX-YY-ZZ   1    0.10          0.0            -3.\nWW-XX-YY-ZZ   1    0.85        180.0            -2.\nWW-XX-YY-ZZ   1    1.35        180.0             1."},
{"N_sp2-C_sp3-O_sb-C_sp3","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.65          0.0             2."},
{"H -C_sp3-O_sb-EP","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"X -C_sp3-O_sb-H ","WW-XX-YY-ZZ   1    0.16          0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1."},
{"X -C_sp3-O_sb-C_sp3","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.1         180.0             2."},
{"X -C_sp3-O_sb-C_sp2","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.80        180.0             1."},
{"X -C_sp3-O_sb-EP","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"C_sp2-C_sp3-O_sb-X ","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.100       180.0             2."},
{"C_sp3-C_sp3-O_sb-X ","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.80        180.0             1."},
{"O_sb-C_sp3-O_sb-X ","WW-XX-YY-ZZ   1    0.10          0.0            -3.\nWW-XX-YY-ZZ   1    0.85        180.0            -2.\nWW-XX-YY-ZZ   1    1.35        180.0             1."},
{"N_sp2-C_sp3-O_sb-X ","WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.65          0.0             2."},
{"H -C_sp3-O_sb-X ","WW-XX-YY-ZZ   1    0.000         0.000           3.0"},
{"X -C_sp3-O_sb-X ","WW-XX-YY-ZZ   3    0.50          0.0             3. "},
{"X -C_sp3-S -X ","WW-XX-YY-ZZ   3    1.00          0.0             3. "},
{"C_sp2-O_sb-P -O_db","WW-XX-YY-ZZ   1    0.0         180.0            -3.\nWW-XX-YY-ZZ   1    0.0         180.0             1."},
{"C_sp2-O_sb-P -O_sb","WW-XX-YY-ZZ   1    0.185181     31.79508        -1.\nWW-XX-YY-ZZ   1    1.256531    351.95960        -2.\nWW-XX-YY-ZZ   1    0.354858    357.24748         3."},
{"C_sp3-O_sb-P -O_sb","WW-XX-YY-ZZ   1    0.25          0.0            -3.\nWW-XX-YY-ZZ   1    1.20          0.0             2."},
{"X -O_sb-P -O_db","WW-XX-YY-ZZ   1    0.0         180.0            -3.\nWW-XX-YY-ZZ   1    0.0         180.0             1."},
{"X -O_sb-P -O_sb","WW-XX-YY-ZZ   1    0.185181     31.79508        -1.\nWW-XX-YY-ZZ   1    1.256531    351.95960        -2.\nWW-XX-YY-ZZ   1    0.354858    357.24748         3."},
{"C_sp2-O_sb-P -X ","WW-XX-YY-ZZ   1    0.185181     31.79508        -1.\nWW-XX-YY-ZZ   1    1.256531    351.95960        -2.\nWW-XX-YY-ZZ   1    0.354858    357.24748         3."},
{"C_sp3-O_sb-P -X ","WW-XX-YY-ZZ   1    0.25          0.0            -3.\nWW-XX-YY-ZZ   1    1.20          0.0             2."},
{"X -O_sb-P -X ","WW-XX-YY-ZZ   3    0.75          0.0             3."},
{"C_sp2-C_sp3-C_sp3-C_sp3","WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1."},
{"H -C_sp3-C_sp3-H ","WW-XX-YY-ZZ   1    0.15          0.0             3. "},
{"H -C_sp3-C_sp3-C_sp3","WW-XX-YY-ZZ   1    0.16          0.0             3."},
{"O_sb-C_sp3-C_sp3-O_sb","WW-XX-YY-ZZ   1    0.144         0.0            -3.\nWW-XX-YY-ZZ   1    1.175         0.0             2."},
{"F -C_sp3-C_sp3-F ","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    1.20        180.0             1."},
{"Cl-C_sp3-C_sp3-Cl","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.45        180.0             1."},
{"Br-C_sp3-C_sp3-Br","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.00        180.0             1."},
{"H -C_sp3-C_sp3-O_sb","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1."},
{"H -C_sp3-C_sp3-F ","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.19          0.0             1."},
{"H -C_sp3-C_sp3-Cl","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1."},
{"H -C_sp3-C_sp3-Br","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.55          0.0             1."},
{"N_sp2-C_sp3-C_sp3-O_sb","WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."},
{"C_sp2-C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1."},
{"H -C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.16          0.0             3."},
{"F -C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    1.20        180.0             1."},
{"Cl-C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.45        180.0             1."},
{"Br-C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.00        180.0             1."},
{"N_sp2-C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."},
{"C_sp3-C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1."},
{"O_sb-C_sp3-C_sp3-X ","WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."},
{"X -C_sp3-C_sp3-C_sp2","WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1."},
{"X -C_sp3-C_sp3-H ","WW-XX-YY-ZZ   1    0.16          0.0             3."},
{"X -C_sp3-C_sp3-F ","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    1.20        180.0             1."},
{"X -C_sp3-C_sp3-Cl","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.45        180.0             1."},
{"X -C_sp3-C_sp3-Br","WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.00        180.0             1."},
{"X -C_sp3-C_sp3-N_sp2","WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."},
{"X -C_sp3-C_sp3-C_sp3","WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1."},
{"X -C_sp3-C_sp3-O_sb","WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."},
{"C_sp2-N_sp2-P -O_db","WW-XX-YY-ZZ   1    0.12          0.0             3."},
{"X -N_sp2-P -O_db","WW-XX-YY-ZZ   1    0.12          0.0             3."},
{"C_sp2-N_sp2-P -X ","WW-XX-YY-ZZ   1    0.12          0.0             3."},
{"X -N_sp2-P -X ","WW-XX-YY-ZZ   1    0.12          0.0             3."},
{"C_sp3-S -S -C_sp3","WW-XX-YY-ZZ   1    3.50          0.0            -2.\nWW-XX-YY-ZZ   1    0.60          0.0             3."},
{"C_sp3-S -S -EP","WW-XX-YY-ZZ   1    0.00          0.0             3."},
{"C_sp3-S -S -X ","WW-XX-YY-ZZ   1    3.50          0.0            -2.\nWW-XX-YY-ZZ   1    0.60          0.0             3."},
{"X -S -S -C_sp3","WW-XX-YY-ZZ   1    3.50          0.0            -2.\nWW-XX-YY-ZZ   1    0.60          0.0             3."},
{"X -S -S -X ","WW-XX-YY-ZZ   1    0.00          0.0             3."},
{"C_sp3-C_sp-C_sp-H ","WW-XX-YY-ZZ   1    0.00          0.0             1."},
{"X -C_sp-C_sp-H ","WW-XX-YY-ZZ   1    0.00          0.0             1."},
{"C_sp3-C_sp-C_sp-X ","WW-XX-YY-ZZ   1    0.00          0.0             1."},
{"X -C_sp-C_sp-X ","WW-XX-YY-ZZ   1    0.00          0.0             1."},
};

std::vector <std::string> DIHEDRAL_SET_DICT_KEYS={
"C_sp2-C_sp2-C_sp2-O_db" ,"C_sp3-C_sp2-C_sp2-C_sp3","C_sp2-C_sp2-C_sp2-X "   ,"O_db-C_sp2-C_sp2-X "    ,"C_sp3-C_sp2-C_sp2-X "   ,"X -C_sp2-C_sp2-X "      ,"H -C_sp2-C_sp3-O_sb"    ,"O_sb-C_sp2-C_sp3-C_sp3",
"C_sp2-C_sp2-C_sp3-H "   ,"N_sp2-C_sp2-C_sp3-N_sp2","N_sp2-C_sp2-C_sp3-C_sp3","O_db-C_sp2-C_sp3-H "    ,"X -C_sp2-C_sp3-O_sb"    ,"X -C_sp2-C_sp3-C_sp3"   ,"X -C_sp2-C_sp3-H "      ,"X -C_sp2-C_sp3-N_sp2"  ,
"H -C_sp2-C_sp3-X "      ,"O_sb-C_sp2-C_sp3-X "    ,"C_sp2-C_sp2-C_sp3-X "   ,"N_sp2-C_sp2-C_sp3-X "   ,"X -C_sp2-C_sp3-X "      ,"O_db-C_sp2-N_sp2-H "    ,"O_db-C_sp2-N_sp2-EP"    ,"C_sp3-C_sp2-N_sp2-EP"  ,
"X -C_sp2-N_sp2-H "      ,"X -C_sp2-N_sp2-EP"      ,"O_db-C_sp2-N_sp2-X "    ,"C_sp3-C_sp2-N_sp2-X "   ,"X -C_sp2-N_sp2-X "      ,"X -C_sp2-O_db-X "       ,"O_db-C_sp2-O_sb-H "     ,"O_db-C_sp2-O_sb-C_sp3" ,
"C_sp2-C_sp2-O_sb-EP"    ,"H -C_sp2-O_sb-EP"       ,"O_db-C_sp2-O_sb-X "     ,"C_sp2-C_sp2-O_sb-X "    ,"H -C_sp2-O_sb-X "       ,"X -C_sp2-O_sb-H "       ,"X -C_sp2-O_sb-C_sp3"    ,"X -C_sp2-O_sb-EP"      ,
"X -C_sp2-O_sb-X "       ,"X -C_sp3-C_sp3-X "      ,"X -C_sp3-C_sp-X "       ,"O_sb-C_sp3-N_sp2-C_sp2" ,"C_sp2-C_sp3-N_sp2-C_sp2","C_sp3-C_sp3-N_sp2-C_sp2","X -C_sp3-N_sp2-C_sp2"   ,"O_sb-C_sp3-N_sp2-X "   ,
"C_sp2-C_sp3-N_sp2-X "   ,"C_sp3-C_sp3-N_sp2-X "   ,"X -C_sp3-N_sp2-X "      ,"C_sp3-C_sp3-N_sp3-C_sp3","H -C_sp3-N_sp3-EP"      ,"C_sp3-C_sp3-N_sp3-EP"   ,"X -C_sp3-N_sp3-C_sp3"   ,"X -C_sp3-N_sp3-EP"     ,
"C_sp3-C_sp3-N_sp3-X "   ,"H -C_sp3-N_sp3-X "      ,"X -C_sp3-N_sp3-X "      ,"C_sp2-C_sp3-O_sb-C_sp3" ,"C_sp3-C_sp3-O_sb-H "    ,"C_sp3-C_sp3-O_sb-C_sp3" ,"C_sp3-C_sp3-O_sb-C_sp2" ,"C_sp3-C_sp3-O_sb-EP"   ,
"O_sb-C_sp3-O_sb-C_sp3"  ,"N_sp2-C_sp3-O_sb-C_sp3" ,"H -C_sp3-O_sb-EP"       ,"X -C_sp3-O_sb-H "       ,"X -C_sp3-O_sb-C_sp3"    ,"X -C_sp3-O_sb-C_sp2"    ,"X -C_sp3-O_sb-EP"       ,"C_sp2-C_sp3-O_sb-X "   ,
"C_sp3-C_sp3-O_sb-X "    ,"O_sb-C_sp3-O_sb-X "     ,"N_sp2-C_sp3-O_sb-X "    ,"H -C_sp3-O_sb-X "       ,"X -C_sp3-O_sb-X "       ,"X -C_sp3-S -X "         ,"C_sp2-O_sb-P -O_db"     ,"C_sp2-O_sb-P -O_sb"    ,
"C_sp3-O_sb-P -O_sb"     ,"X -O_sb-P -O_db"        ,"X -O_sb-P -O_sb"        ,"C_sp2-O_sb-P -X "       ,"C_sp3-O_sb-P -X "       ,"X -O_sb-P -X "          ,"C_sp2-C_sp3-C_sp3-C_sp3","H -C_sp3-C_sp3-H "     ,
"H -C_sp3-C_sp3-C_sp3"   ,"O_sb-C_sp3-C_sp3-O_sb"  ,"F -C_sp3-C_sp3-F "      ,"Cl-C_sp3-C_sp3-Cl"      ,"Br-C_sp3-C_sp3-Br"      ,"H -C_sp3-C_sp3-O_sb"    ,"H -C_sp3-C_sp3-F "      ,"H -C_sp3-C_sp3-Cl"     ,
"H -C_sp3-C_sp3-Br"      ,"N_sp2-C_sp3-C_sp3-O_sb" ,"C_sp2-C_sp3-C_sp3-X "   ,"H -C_sp3-C_sp3-X "      ,"F -C_sp3-C_sp3-X "      ,"Cl-C_sp3-C_sp3-X "      ,"Br-C_sp3-C_sp3-X "      ,"N_sp2-C_sp3-C_sp3-X "  ,
"C_sp3-C_sp3-C_sp3-X "   ,"O_sb-C_sp3-C_sp3-X "    ,"X -C_sp3-C_sp3-C_sp2"   ,"X -C_sp3-C_sp3-H "      ,"X -C_sp3-C_sp3-F "      ,"X -C_sp3-C_sp3-Cl"      ,"X -C_sp3-C_sp3-Br"      ,"X -C_sp3-C_sp3-N_sp2"  ,
"X -C_sp3-C_sp3-C_sp3"   ,"X -C_sp3-C_sp3-O_sb"    ,"C_sp2-N_sp2-P -O_db"    ,"X -N_sp2-P -O_db"       ,"C_sp2-N_sp2-P -X "      ,"X -N_sp2-P -X "         ,"C_sp3-S -S -C_sp3"      ,"C_sp3-S -S -EP"        ,
"C_sp3-S -S -X "         ,"X -S -S -C_sp3"         ,"X -S -S -X "            ,"C_sp3-C_sp-C_sp-H "     ,"X -C_sp-C_sp-H "        ,"C_sp3-C_sp-C_sp-X "     ,"X -C_sp-C_sp-X "        ,"X -C_sp2-S-X "         ,
};


std::map <std::string,std::string> dihedral_library={
{"X -C -C -X ","X -C -C -X     4   14.500       180.000           2.000"},
{"X -C -C*-X ","X -C -C*-X     4   14.500       180.000           2.000"},
{"X -C -CA-X ","X -C -CA-X     4   14.500       180.000           2.000"},
{"X -C -CB-X ","X -C -CB-X     4   14.500       180.000           2.000"},
{"X -C -CD-X ","X -C -CD-X     4   14.500       180.000           2.000"},
{"X -C -CK-X ","X -C -CK-X     4   14.500       180.000           2.000"},
{"X -C -CM-X ","X -C -CM-X     4   14.500       180.000           2.000"},
{"X -C -CN-X ","X -C -CN-X     4   14.500       180.000           2.000"},
{"X -C -CQ-X ","X -C -CQ-X     4   14.500       180.000           2.000"},
{"X -C -CT-X ","X -C -CT-X     1    0.000         0.000          -4.000\nX -C -CT-X     1    0.420         0.000          -3.000\nX -C -CT-X     1    0.270         0.000          -2.000\nX -C -CT-X     1    0.000         0.000           1.000"},
{"X -C -CZ-X ","X -C -CZ-X     4   14.500       180.000           2.000"},
{"X -C -H1-X ","X -C -H1-X     1    0.000         0.000          -4.000\nX -C -H1-X     1    0.420         0.000          -3.000\nX -C -H1-X     1    0.270         0.000          -2.000\nX -C -H1-X     1    0.000         0.000           1.000"},
{"X -C -HC-X ","X -C -HC-X     1    0.000         0.000          -4.000\nX -C -HC-X     1    0.420         0.000          -3.000\nX -C -HC-X     1    0.270         0.000          -2.000\nX -C -HC-X     1    0.000         0.000           1.000"},
{"X -C -N -X ","X -C -N -X     4   14.500       180.000           2.000"},
{"X -C -N2-X ","X -C -N2-X     4   14.500       180.000           2.000"},
{"X -C -NB-X ","X -C -NB-X     4    5.400       180.000           2.000\nX -C -NB-X     4   14.500       180.000           2.000"},
{"X -C -NT-X ","X -C -NT-X     1    0.000         0.000          -4.000\nX -C -NT-X     1    0.420         0.000          -3.000\nX -C -NT-X     1    0.270         0.000          -2.000\nX -C -NT-X     1    0.000         0.000           1.000"},
{"X -C*-C*-X ","X -C*-C*-X     4   14.500       180.000           2.000"},
{"X -C*-CA-X ","X -C*-CA-X     4   14.500       180.000           2.000"},
{"X -C*-CC-X ","X -C*-CC-X     4   14.500       180.000           2.000"},
{"X -C*-CM-X ","X -C*-CM-X     4   14.500       180.000           2.000"},
{"X -C*-CQ-X ","X -C*-CQ-X     4   14.500       180.000           2.000"},
{"X -CA-CA-X ","X -CA-CA-X     4   14.500       180.000           2.000"},
{"X -CA-CB-X ","X -CA-CB-X     4    6.600       180.000           2.000"},
{"X -CA-CC-X ","X -CA-CC-X     4   14.500       180.000           2.000"},
{"X -CA-CD-X ","X -CA-CD-X     4   14.500       180.000           2.000"},
{"X -CA-CK-X ","X -CA-CK-X     4   14.500       180.000           2.000"},
{"X -CA-CM-X ","X -CA-CM-X     4   14.500       180.000           2.000"},
{"X -CA-CQ-X ","X -CA-CQ-X     4   14.500       180.000           2.000"},
{"X -CA-CZ-X ","X -CA-CZ-X     4   14.500       180.000           2.000"},
{"X -CA-N -X ","X -CA-N -X     4   14.500       180.000           2.000"},
{"X -CA-N*-X ","X -CA-N*-X     4   14.500       180.000           2.000"},
{"X -CA-N2-X ","X -CA-N2-X     4   14.500       180.000           2.000"},
{"X -CA-NB-X ","X -CA-NB-X     4    6.000       180.000           2.000\nX -CA-NB-X     4   14.500       180.000           2.000"},
{"X -CA-NC-X ","X -CA-NC-X     4   14.500       180.000           2.000"},
{"X -CA-OS-X ","X -CA-OS-X     4   14.500       180.000           2.000"},
{"X -CB-CM-X ","X -CB-CM-X     4   14.500       180.000           2.000"},
{"X -CB-CT-X ","X -CB-CT-X     1    0.000         0.000          -4.000\nX -CB-CT-X     1    0.420         0.000          -3.000\nX -CB-CT-X     1    0.270         0.000          -2.000\nX -CB-CT-X     1    0.000         0.000           1.000"},
{"X -CB-N -X ","X -CB-N -X     4    6.600       180.000           2.000"},
{"X -CB-N*-X ","X -CB-N*-X     4   14.500       180.000           2.000"},
{"X -CB-NB-X ","X -CB-NB-X     4   14.500       180.000           2.000"},
{"X -CB-OS-X ","X -CB-OS-X     1    0.000         0.000          -4.000\nX -CB-OS-X     1    0.420         0.000          -3.000\nX -CB-OS-X     1    0.270         0.000          -2.000\nX -CB-OS-X     1    0.000         0.000           1.000"},
{"X -CC-CM-X ","X -CC-CM-X     4   14.500       180.000           2.000"},
{"X -CC-N -X ","X -CC-N -X     4   14.500       180.000           2.000"},
{"X -CD-CD-X ","X -CD-CD-X     4   14.500       180.000           2.000"},
{"X -CD-CM-X ","X -CD-CM-X     4   14.500       180.000           2.000"},
{"X -CD-CQ-X ","X -CD-CQ-X     4   14.500       180.000           2.000"},
{"X -CD-N -X ","X -CD-N -X     4   14.500       180.000           2.000"},
{"X -CD-NB-X ","X -CD-NB-X     4   14.500       180.000           2.000"},
{"X -CK-CM-X ","X -CK-CM-X     4   14.500       180.000           2.000"},
{"X -CK-CT-X ","X -CK-CT-X     6    0.000         0.000           2.000"},
{"X -CK-N -X ","X -CK-N -X     4   14.500       180.000           2.000"},
{"X -CK-N*-X ","X -CK-N*-X     4   14.500       180.000           2.000"},
{"X -CK-NB-X ","X -CK-NB-X     4   14.500       180.000           2.000"},
{"X -CM-CM-X ","X -CM-CM-X     4   14.500       180.000           2.000"},
{"X -CM-CQ-X ","X -CM-CQ-X     4   14.500       180.000           2.000"},
{"X -CM-CT-X ","X -CM-CT-X     1    0.000         0.000          -4.000\nX -CM-CT-X     1    0.420         0.000          -3.000\nX -CM-CT-X     1    0.270         0.000          -2.000\nX -CM-CT-X     1    0.000         0.000           1.000"},
{"X -CM-CZ-X ","X -CM-CZ-X     4   14.500       180.000           2.000"},
{"X -CM-N -X ","X -CM-N -X     4   14.500       180.000           2.000"},
{"X -CM-N*-X ","X -CM-N*-X     4   14.500       180.000           2.000"},
{"X -CM-N2-X ","X -CM-N2-X     4   14.500       180.000           2.000"},
{"X -CM-NA-X ","X -CM-NA-X     2    8.000       180.000           2.000"},
{"X -CM-NB-X ","X -CM-NB-X     4    7.400       180.000           2.000\nX -CM-NB-X     4   14.500       180.000           2.000"},
{"X -CM-NC-X ","X -CM-NC-X     4   14.500       180.000           2.000"},
{"X -CM-NT-X ","X -CM-NT-X     1    0.000         0.000          -4.000\nX -CM-NT-X     1    0.420         0.000          -3.000\nX -CM-NT-X     1    0.270         0.000          -2.000\nX -CM-NT-X     1    0.000         0.000           1.000"},
{"X -CM-OH-X ","X -CM-OH-X     1    0.000         0.000          -4.000\nX -CM-OH-X     1    0.420         0.000          -3.000\nX -CM-OH-X     1    0.270         0.000          -2.000\nX -CM-OH-X     1    0.000         0.000           1.000"},
{"X -CM-OS-X ","X -CM-OS-X     1    0.000         0.000          -4.000\nX -CM-OS-X     1    0.420         0.000          -3.000\nX -CM-OS-X     1    0.270         0.000          -2.000\nX -CM-OS-X     1    0.000         0.000           1.000"},
{"X -CQ-CT-X ","X -CQ-CT-X     1    0.000         0.000          -4.000\nX -CQ-CT-X     1    0.420         0.000          -3.000\nX -CQ-CT-X     1    0.270         0.000          -2.000\nX -CQ-CT-X     1    0.000         0.000           1.000"},
{"X -CQ-N -X ","X -CQ-N -X     4   14.500       180.000           2.000"},
{"X -CQ-N2-X ","X -CQ-N2-X     4   14.500       180.000           2.000"},
{"X -CQ-NB-X ","X -CQ-NB-X     4   14.500       180.000           2.000"},
{"X -CQ-NC-X ","X -CQ-NC-X     4   14.500       180.000           2.000"},
{"X -CQ-OS-X ","X -CQ-OS-X     1    0.000         0.000          -4.000\nX -CQ-OS-X     1    0.420         0.000          -3.000\nX -CQ-OS-X     1    0.270         0.000          -2.000\nX -CQ-OS-X     1    0.000         0.000           1.000"},
{"X -CR-N -X ","X -CR-N -X     4    9.300       180.000           2.000"},
{"X -CR-N*-X ","X -CR-N*-X     4   14.500       180.000           2.000"},
{"X -CR-N2-X ","X -CR-N2-X     4   14.500       180.000           2.000"},
{"X -CT-N -X ","X -CT-N -X     1    0.000         0.000          -4.000\nX -CT-N -X     1    0.420         0.000          -3.000\nX -CT-N -X     1    0.270         0.000          -2.000\nX -CT-N -X     1    0.000         0.000           1.000"},
{"X -CT-N2-X ","X -CT-N2-X     1    0.000         0.000          -4.000\nX -CT-N2-X     1    0.420         0.000          -3.000\nX -CT-N2-X     1    0.270         0.000          -2.000\nX -CT-N2-X     1    0.000         0.000           1.000"},
{"X -CT-NB-X ","X -CT-NB-X     6    0.000         0.000           2.000"},
{"X -CT-P -X ","X -CT-P -X     1    0.000         0.000          -4.000\nX -CT-P -X     1    0.420         0.000          -3.000\nX -CT-P -X     1    0.270         0.000          -2.000\nX -CT-P -X     1    0.000         0.000           1.000"},
{"X -CZ-N -X ","X -CZ-N -X     4   14.500       180.000           2.000"},
{"X -CZ-NB-X ","X -CZ-NB-X     4   14.500       180.000           2.000"},
{"X -H -NB-X ","X -H -NB-X     1    0.000         0.000          -4.000\nX -H -NB-X     1    0.420         0.000          -3.000\nX -H -NB-X     1    0.270         0.000          -2.000\nX -H -NB-X     1    0.000         0.000           1.000"},
{"X -HP-NB-X ","X -HP-NB-X     1    0.000         0.000          -4.000\nX -HP-NB-X     1    0.420         0.000          -3.000\nX -HP-NB-X     1    0.270         0.000          -2.000\nX -HP-NB-X     1    0.000         0.000           1.000"},
{"X -N -OS-X ","X -N -OS-X     1    0.000         0.000          -4.000\nX -N -OS-X     1    0.420         0.000          -3.000\nX -N -OS-X     1    0.270         0.000          -2.000\nX -N -OS-X     1    0.000         0.000           1.000"},
{"X -N*-NB-X ","X -N*-NB-X     4   14.500       180.000           2.000"},
{"X -N2-NB-X ","X -N2-NB-X     4   14.500       180.000           2.000"},
{"X -N3-NB-X ","X -N3-NB-X     1    0.000         0.000          -4.000\nX -N3-NB-X     1    0.420         0.000          -3.000\nX -N3-NB-X     1    0.270         0.000          -2.000\nX -N3-NB-X     1    0.000         0.000           1.000"},
{"X -NB-NT-X ","X -NB-NT-X     1    0.000         0.000          -4.000\nX -NB-NT-X     1    0.420         0.000          -3.000\nX -NB-NT-X     1    0.270         0.000          -2.000\nX -NB-NT-X     1    0.000         0.000           1.000"},
{"X -NB-O -X ","X -NB-O -X     1    0.000         0.000          -4.000\nX -NB-O -X     1    0.420         0.000          -3.000\nX -NB-O -X     1    0.270         0.000          -2.000\nX -NB-O -X     1    0.000         0.000           1.000"},
{"X -OS-SO-X ","X -OS-SO-X     1    0.000         0.000          -4.000\nX -OS-SO-X     1    0.420         0.000          -3.000\nX -OS-SO-X     1    0.270         0.000          -2.000\nX -OS-SO-X     1    0.000         0.000           1.000"},
{"C -C*-S -CB","C -C*-S -CB   2    4.00        180.000           2. "},
{"C*-C*-S -CB","C*-C*-S -CB   2    4.00        180.000           2. "},
{"X -C*-S -X ","X -C*-S -X    2    4.00        180.000           2. "},
{"C*-C*-N*-CT","C*-C*-N*-CT   2    7.40        180.000           2. "},
{"C*-C*-N*-CB","C*-C*-N*-CB   2    7.40        180.000           2. "},
{"C -C*-N*-CB","C -C*-N*-CB   2    7.40        180.000           2. "},
{"C -C*-N*-CT","C -C*-N*-CT   2    7.40        180.000           2. "},
{"X -C*-N*-X ","X -C*-N*-X    2    7.40        180.000           2. "},
{"CA-CB-S -C*","CA-CB-S -C*   2    4.00        180.000           2. "},
{"CB-CB-S -C*","CB-CB-S -C*   2    4.00        180.000           2. "},
{"X -CB-S -X ","X -CB-S -X    2    4.00        180.000           2. "},
};

std::vector <std::string> dihedral_library_keys={
"X -C -C -X ","X -C -C*-X ","X -C -CA-X ","X -C -CB-X ","X -C -CD-X ","X -C -CK-X ","X -C -CM-X ","X -C -CN-X ","X -C -CQ-X ","X -C -CT-X ",
"X -C -CZ-X ","X -C -H1-X ","X -C -HC-X ","X -C -N -X ","X -C -N2-X ","X -C -NB-X ","X -C -NT-X ","X -C*-C*-X ","X -C*-CA-X ","X -C*-CC-X ",
"X -C*-CM-X ","X -C*-CQ-X ","X -CA-CA-X ","X -CA-CB-X ","X -CA-CC-X ","X -CA-CD-X ","X -CA-CK-X ","X -CA-CM-X ","X -CA-CQ-X ","X -CA-CZ-X ",
"X -CA-N -X ","X -CA-N*-X ","X -CA-N2-X ","X -CA-NB-X ","X -CA-NC-X ","X -CA-OS-X ","X -CB-CM-X ","X -CB-CT-X ","X -CB-N -X ","X -CB-N*-X ",
"X -CB-NB-X ","X -CB-OS-X ","X -CC-CM-X ","X -CC-N -X ","X -CD-CD-X ","X -CD-CM-X ","X -CD-CQ-X ","X -CD-N -X ","X -CD-NB-X ","X -CK-CM-X ",
"X -CK-CT-X ","X -CK-N -X ","X -CK-N*-X ","X -CK-NB-X ","X -CM-CM-X ","X -CM-CQ-X ","X -CM-CT-X ","X -CM-CZ-X ","X -CM-N -X ","X -CM-N*-X ",
"X -CM-N2-X ","X -CM-NA-X ","X -CM-NB-X ","X -CM-NC-X ","X -CM-NT-X ","X -CM-OH-X ","X -CM-OS-X ","X -CQ-CT-X ","X -CQ-N -X ","X -CQ-N2-X ",
"X -CQ-NB-X ","X -CQ-NC-X ","X -CQ-OS-X ","X -CR-N -X ","X -CR-N*-X ","X -CR-N2-X ","X -CT-N -X ","X -CT-N2-X ","X -CT-NB-X ","X -CT-P -X ",
"X -CZ-N -X ","X -CZ-NB-X ","X -H -NB-X ","X -HP-NB-X ","X -N -OS-X ","X -N*-NB-X ","X -N2-NB-X ","X -N3-NB-X ","X -NB-NT-X ","X -NB-O -X ",
"X -OS-SO-X ","C -C*-S -CB","C*-C*-S -CB","X -C*-S -X ","C*-C*-N*-CT","C*-C*-N*-CB","C -C*-N*-CB","C -C*-N*-CT","X -C*-N*-X ","CA-CB-S -C*",
"CB-CB-S -C*","X -CB-S -X ",
};

// vdW library (map of bond-pair strings to full frcmod string)
std::map <std::string,std::string> vdw_library={
{"H ","H           0.6000  0.0157  "},{"HO","HO          0.0000  0.0000  "},{"HS","HS          0.6000  0.0157  "},{"HC","HC          1.4870  0.0157  "},
{"H1","H1          1.3870  0.0157  "},{"H2","H2          1.2870  0.0157  "},{"H3","H3          1.1870  0.0157  "},{"HP","HP          1.1000  0.0157  "},
{"HA","HA          1.4590  0.0150  "},{"H4","H4          1.4090  0.0150  "},{"H5","H5          1.3590  0.0150  "},{"HW","HW          0.0000  0.0000  "},
{"HZ","HZ          1.4590  0.0150  "},{"O ","O           1.6612  0.2100  "},{"O2","O2          1.6612  0.2100  "},{"OW","OW          1.7683  0.1520  "},
{"OH","OH          1.7210  0.2104  "},{"OS","OS          1.6837  0.1700  "},{"OP","OP          1.8500  0.1700  "},{"C*","C*          1.9080  0.0860  "},
{"CI","CI          1.9080  0.1094  "},{"C5","C5          1.9080  0.0860  "},{"C4","C4          1.9080  0.0860  "},{"CT","CT          1.9080  0.1094  "},{"CG","CG          1.9080  0.1094  "},
{"CX","CX          1.9080  0.1094  "},{"C ","C           1.9080  0.0860  "},{"CA","CA          1.9080  0.0860  "},{"CB","CB          1.9080  0.0860  "},
{"CC","CC          1.9080  0.0860  "},{"CK","CK          1.9080  0.0860  "},{"CM","CM          1.9080  0.0860  "},{"CQ","CQ          1.9080  0.0860  "},
{"CR","CR          1.9080  0.0860  "},{"CZ","CZ          1.9080  0.0860  "},{"N ","N           1.8240  0.1700  "},{"N3","N3          1.8240  0.1700  "},
{"N*","N*          1.8240  0.1700  "},{"N2","N2          1.8240  0.1700  "},{"NA","NA          1.8240  0.1700  "},{"NB","NB          1.8240  0.1700  "},
{"NC","NC          1.8240  0.1700  "},{"NT","NT          1.8240  0.1700  "},{"NY","NY          1.8240  0.1700  "},{"S ","S           2.0000  0.2500  "},
{"SG","SG          2.0000  0.2500  "},{"SH","SH          2.0000  0.2500  "},{"SO","SO          2.0000  0.2500  "},{"SS","SS          2.0000  0.2500  "},
{"P ","P           2.1000  0.2000  "},{"MG","MG          0.7926  0.8947  "},{"C0","C0          1.7131  0.459789"},{"F ","F           1.75    0.061   "},
{"Cl","Cl          1.948   0.265   "},{"Br","Br          2.22    0.320   "},{"I ","I           2.35    0.40    "},{"N1","N1          1.8240  0.1700  "}};


std::vector <std::string> vdw_library_keys={
{"H "},{"HO"},{"HS"},{"HC"},{"H1"},{"H2"},{"H3"},{"HP"},{"HA"},{"H4"},
{"H5"},{"HW"},{"HZ"},{"O "},{"O2"},{"OW"},{"OH"},{"OS"},{"OP"},{"C*"},
{"CI"},{"C5"},{"C4"},{"CT"},{"CG"},{"CX"},{"C "},{"CA"},{"CB"},{"CC"},{"CK"},
{"CM"},{"CQ"},{"CR"},{"CZ"},{"N "},{"N3"},{"N*"},{"N2"},{"NA"},{"NB"},
{"NC"},{"NT"},{"NY"},{"S "},{"SG"},{"SH"},{"SO"},{"SS"},{"P "},{"MG"},
{"C0"},{"F "},{"Cl"},{"Br"},{"I "},{"N1"},};

std::string FindMissingMass(std::string atom_name)
{
    for (int i=0; i < mass_library_keys.size(); i++)
    {
        if (mass_library_keys[i] == atom_name)
        {
            std::string result = mass_library[atom_name] + "\n";
            // std::cout << result << std::endl;
            return result;
        }
    }
    // if mass is not found, return NewMass with value of "MISSING INFORMATION"
    std::stringstream buffer;
    buffer.str("");
    buffer << std::setw(2) << std::left << atom_name << " MISSING INFORMATION\n";
    return buffer.str();
}

std::string FindMissingBond(std::string a1, std::string a2)
{
    std::stringstream buffer;
    std::string pair1;
    std::string pair2;
    //check if 'a1-a2' exists in the bond library.
    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    pair1 = buffer.str();

    //check if 'a2-a1' exists in the bond library.
    buffer.str("");
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a1;
    pair2 = buffer.str();

    for (int i=0; i < bond_library_keys.size(); i++)
    {
        if (bond_library_keys[i] == pair1)
        {
            return bond_library[pair1];
        }
        if (bond_library_keys[i] == pair2)
        {
            return bond_library[pair2];
        }
    }

    // If neither version of the bond listing are present, obtain basic atom types and generate blurry string
    std::string base1 = BASIC_ATOM_TYPES[a1];
    std::string base2 = BASIC_ATOM_TYPES[a2];

    pair1 = base1+"-"+base2;
    pair2 = base2+"-"+base1;
   
    for (int i=0; i < bond_defaults_keys.size(); i++)
    {
        if (bond_defaults_keys[i] == pair1)
        {
            buffer.str("");
            buffer << std::setw(2) << std::left << a1;
            buffer << "-";
            buffer << std::setw(2) << std::left << a2;
            buffer << bond_defaults[pair1] << std::endl;
            return buffer.str();
        }
        if (bond_defaults_keys[i] == pair2)
        {
            buffer.str("");
            buffer << std::setw(2) << std::left << a2;
            buffer << "-";
            buffer << std::setw(2) << std::left << a1;
            buffer << bond_defaults[pair2] << std::endl;
            return buffer.str();
        }
    }
    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << " MISSING INFORMATION " << std::endl;
    return buffer.str();
}

std::string FindMissingAngle(std::string a1, std::string a2, std::string a3)
{
    std::stringstream buffer;
    std::string angle1;
    std::string angle2;
    //check if 'a1-a2-a3' exists in the bond library.
    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a3;
    angle1 = buffer.str();

    //check if 'a3-a2-a1' exists in the bond library.
    buffer.str("");
    buffer << std::setw(2) << std::left << a3;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a1;
    angle2 = buffer.str();

    for (int i=0; i < angle_library_keys.size(); i++)
    {
        if (angle_library_keys[i] == angle1)
        {
            return angle_library[angle1];
        }
        if (angle_library_keys[i] == angle2)
        {
            return angle_library[angle2];
        }
    }
    // No exact match found, expand outward
    // If neither version of the bond listing are present, obtain basic atom types and generate blurry string
    // std::cout << "Unable to find exact match to " << angle1 << " or " << angle2 << std::endl;
    std::string base1 = BASIC_ATOM_TYPES[a1];
    std::string base2 = BASIC_ATOM_TYPES[a2];
    std::string base3 = BASIC_ATOM_TYPES[a3];

    angle1 = base1 + "-" + base2 + "-" + base3;
    angle2 = base3 + "-" + base2 + "-" + base1;
    // std::cout << "searching for match to " << angle1 << " or " << angle2 << std::endl;
   
    for (int i=0; i < angle_defaults_keys.size(); i++)
    {
        if (angle_defaults_keys[i] == angle1)
        {
            buffer.str("");
            buffer << std::setw(2) << std::left << a1;
            buffer << "-";
            buffer << std::setw(2) << std::left << a2;
            buffer << "-";
            buffer << std::setw(2) << std::left << a3;
            buffer << angle_defaults[angle1] << std::endl;
            return buffer.str();
        }
        if (angle_defaults_keys[i] == angle2)
        {
            buffer.str("");
            buffer << std::setw(2) << std::left << a3;
            buffer << "-";
            buffer << std::setw(2) << std::left << a2;
            buffer << "-";
            buffer << std::setw(2) << std::left << a1;
            buffer << angle_defaults[angle2] << std::endl;
            return buffer.str();
        }
    }
    angle1 = base2;
    for (int i=0; i < angle_defaults_keys.size(); i++)
    {
        if (angle_defaults_keys[i] == angle1)
        {
            buffer.str("");
            buffer << std::setw(2) << std::left << a1;
            buffer << "-";
            buffer << std::setw(2) << std::left << a2;
            buffer << "-";
            buffer << std::setw(2) << std::left << a3;
            buffer << angle_defaults[angle1] << std::endl;
            return buffer.str();
        }
    }
    buffer << " MISSING INFORMATION " << std::endl;
    return buffer.str();
    
}

std::string FindMissingDihedral(std::string a1, std::string a2, std::string a3, std::string a4)
{
    std::stringstream buffer;
    std::string dihedral1;
    std::string dihedral2;
    //check if 'a1-a2-a3-a4' exists in the bond library.
    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a3;
    buffer << "-";
    buffer << std::setw(2) << std::left << a4;
    dihedral1 = buffer.str();

    //check if 'a4-a3-a2-a1' exists in the bond library.
    buffer.str("");
    buffer << std::setw(2) << std::left << a4;
    buffer << "-";
    buffer << std::setw(2) << std::left << a3;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a1;
    dihedral2 = buffer.str();

    for (int i=0; i < dihedral_library_keys.size(); i++)
    {
        if (dihedral_library_keys[i] == dihedral1)
        {
            return dihedral_library[dihedral1];
        }
        if (dihedral_library_keys[i] == dihedral2)
        {
            return dihedral_library[dihedral2];
        }
    }
    std::string replace1 = dihedral1;
    std::string replace2 = dihedral2;

    // No exact match found, expand outward
    // If neither version of the dihedral listing are present, 
    // obtain basic atom types and generate blurry strings
    // std::cout << "Unable to find exact match to " << dihedral1 << " or " << dihedral2 << std::endl;
    std::string base1 = BASIC_ATOM_TYPES[a1];
    std::string base2 = BASIC_ATOM_TYPES[a2];
    std::string base3 = BASIC_ATOM_TYPES[a3];
    std::string base4 = BASIC_ATOM_TYPES[a4];

    dihedral1 = "X -" + a2 + "-" + a3 + "-X ";
    dihedral2 = "X -" + a3 + "-" + a2 + "-X ";
    // std::cout << "searching for match to " << dihedral1 << " or " << dihedral2 << std::endl;
    for (int i=0; i < dihedral_library_keys.size(); i++)
    {
        if (dihedral_library_keys[i] == dihedral1)
        {
            return dihedral_library[dihedral1];
        }
        if (dihedral_library_keys[i] == dihedral2)
        {
            return dihedral_library[dihedral2];
        }
    }


    // substring_replace(DIHEDRAL_SET_DICT[dihedral1], "WW-XX-YY-ZZ", replace1)
    // substring_replace(DIHEDRAL_SET_DICT[dihedral2], "WW-XX-YY-ZZ", replace2)

    dihedral1 = base1 + "-" + base2 + "-" + base3 + "-" + base4;
    dihedral2 = base4 + "-" + base3 + "-" + base2 + "-" + base1;
    // std::cout << "searching for match to " << dihedral1 << " or " << dihedral2 << std::endl;
    for (int i=0; i < DIHEDRAL_SET_DICT_KEYS.size(); i++)
    {
        if (DIHEDRAL_SET_DICT_KEYS[i] == dihedral1)
        {
            return substring_replace(DIHEDRAL_SET_DICT[dihedral1], "WW-XX-YY-ZZ", replace1);
        }
        if (DIHEDRAL_SET_DICT_KEYS[i] == dihedral2)
        {
            return substring_replace(DIHEDRAL_SET_DICT[dihedral2], "WW-XX-YY-ZZ", replace2);
        }
    }

    dihedral1 = "X -" + base2 + "-" + base3 + "-" + base4;
    dihedral2 = "X -" + base3 + "-" + base2 + "-" + base1;
    // std::cout << "searching for match to " << dihedral1 << " or " << dihedral2 << std::endl;
    for (int i=0; i < DIHEDRAL_SET_DICT_KEYS.size(); i++)
    {
        if (DIHEDRAL_SET_DICT_KEYS[i] == dihedral1)
        {
            return substring_replace(DIHEDRAL_SET_DICT[dihedral1], "WW-XX-YY-ZZ", replace1);
        }
        if (DIHEDRAL_SET_DICT_KEYS[i] == dihedral2)
        {
            return substring_replace(DIHEDRAL_SET_DICT[dihedral2], "WW-XX-YY-ZZ", replace2);
        }
    }

    dihedral1 = base1 + "-" + base2 + "-" + base3 + "-X ";
    dihedral2 = base4 + "-" + base3 + "-" + base2 + "-X ";
    // std::cout << "searching for match to " << dihedral1 << " or " << dihedral2 << std::endl;
    for (int i=0; i < DIHEDRAL_SET_DICT_KEYS.size(); i++)
    {
        if (DIHEDRAL_SET_DICT_KEYS[i] == dihedral1)
        {
            return substring_replace(DIHEDRAL_SET_DICT[dihedral1], "WW-XX-YY-ZZ", replace1);
        }
        if (DIHEDRAL_SET_DICT_KEYS[i] == dihedral2)
        {
            return substring_replace(DIHEDRAL_SET_DICT[dihedral2], "WW-XX-YY-ZZ", replace2);
        }
    }

    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a3;
    buffer << "-";
    buffer << std::setw(2) << std::left << a4;
    buffer << " MISSING INFORMATION " << std::endl;
    return buffer.str();
}

std::string FindMissingTorsion(std::string a1, std::string a2, std::string a3, std::string a4)
{
    std::stringstream buffer;
    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << "-";
    buffer << std::setw(2) << std::left << a2;
    buffer << "-";
    buffer << std::setw(2) << std::left << a3;
    buffer << "-";
    buffer << std::setw(2) << std::left << a4;
    buffer << " MISSING INFORMATION " << std::endl;
    return buffer.str();
    
}

std::string FindMissingNonbonded(std::string a1)
{

    for (int i=0; i < vdw_library_keys.size(); i++)
    {
        if (vdw_library_keys[i] == a1)
        {
            return vdw_library[a1];
        }
    }

    std::stringstream buffer;
    buffer.str("");
    buffer << std::setw(2) << std::left << a1;
    buffer << " MISSING INFORMATION " << std::endl;
    return buffer.str();
}
