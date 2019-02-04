//Gtest Includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

//Filter Include
#include "port_despike.h"

//Other Includes
#include <iostream>

namespace{

using namespace Eigen;
using namespace despike;

class DespikeTest:public ::testing::Test
{
public:
    DespikeTest(): generator(std::random_device{}()) {};
    ~DespikeTest() override {}
    void TearDown() override {}
    void SetUp() override {}

    std::mt19937 generator;

    //Define some basic input parameters for all tests
    int nSamples = 1001; //number of samples in time series.
    int nTerms=32;
    int nCoef=2*nTerms+1;

    //Filter parameters
    double samplerate = 64.; //Samplerate
    double nyquist = samplerate/2; //Nyquist Sampling rate

    int nSpikes = 10;
    int SpikeSigma = 30;
    int SpikeLim = 5;
    int despikeWindow = nCoef;

    double timeConstant = 0.015;   //assumed time constant [seconds]

    bool isLowpassed = 0;
    bool isDespiked = 0;

    //Error tolerance for google test
    double err = 1e-5;

    //Parameters for random time series
    const double distmean = 0.0;
    const double diststddev = 0.01;

    //std::vector<float> no_spike = {-0.00146382, 0.0013453, -0.0187138, 0.0046065, -0.00214253, 0.00163712, -0.00827944, 0.00298595, 0.0105547, 0.000102154, 0.0117457, -0.00546841, -0.0104944, 0.00660682, -0.00625276, 0.0148596, -0.00829081, -0.0255912, -0.00888707, -0.00539781, 0.0101922, -0.00628956, -0.00482589, 0.00339587, -0.00121306, 0.0210886, -0.00371003, -0.00287389, -0.0230144, -0.0105935, -0.000615274, 0.0145502, 0.0135433, 0.00925328, -0.00243275, 0.0151561, 0.00197497, 0.0100886, 0.00439499, 0.00438945, 0.00645743, -0.00128149, -0.0168599, 0.0177643, -0.00613857, 0.00469861, -0.00582398, 0.00668493, -0.00103692, 0.00149386, 0.00624049, 0.0153727, 0.0117067, 0.0107825, -0.0205006, 0.0117196, -0.0145473, 0.00136395, -0.0111552, -0.0171463, 0.0112422, -0.0173985, -0.0147975, -0.0158694, 0.0148247, -0.00727862, 0.00754843, -0.001128, 0.00984235, 0.00326633, -0.0103745, -0.000764704, -0.0208402, 0.00389231, 0.00243215, 0.00455092, 0.00275194, 0.0291628, 0.00272422, -0.0320464, 0.0186225, -0.0209501, 0.0105544, 0.00310367, -1.22802e-05, 0.00404831, -0.0108115, 0.0141863, -0.00400148, 0.00926096, -0.00358203, 0.00126072, 0.00387892, -0.00569566, -0.00634654, 0.00882249, -0.00677104, 0.00204175, 0.0135715, -0.02453, -0.00315325, -0.00379922, -0.00608541, 0.0135717, 0.000195746, 0.0132359, -0.000912438, -0.00208138, -0.0161209, 0.00281664, 0.00785215, -0.00316253, 0.00353801, -0.00271609, -0.0177443, -0.000590157, -0.0153723, -0.000539041, 0.00386642, 0.00129153, 0.0154308, -0.00388037, 0.000763218, -0.00644921, -0.0168173, -0.00900313, -0.00383684, -0.00587587, -0.00130386, -0.0237398, 0.0101014, -0.00981124, 0.00465062, 0.0117503, -0.0154914, -0.00944311, -0.00524505, 0.0106121, 0.0152626, -0.0108715, -0.00395067, -0.018694, 0.00886654, -0.000245811, 0.00281015, -0.0043059, -0.00282224, 0.0098182, -0.00226851, 0.00680064, 0.00137044, -0.011258, -0.00171483, -0.000395065, 0.0151453, -0.00126506, 0.000306945, -0.00368253, 0.00983281, -0.00149726, -0.00584821, -0.00432586, 0.0141966, -0.00127846, -0.00704493, 0.0024827, 0.0181503, 0.01151, 0.0132445, -0.0084942, 0.0101566, 0.0116719, -0.000821409, -0.000890962, -0.000626403, -0.00568668, 0.00450628, -0.00503789, 0.00950093, -0.00502095, 0.0182362, -0.0093429, 0.00246156, -0.0033775, 0.00884862, -0.00787764, -0.00248331, -0.0038751, -0.00916308, -0.0128537, 0.00527725, -0.013623, -0.0124964, -0.00396654, -0.00314504, -0.0119639, -0.0140579, 0.00954706, -0.00564152, 0.00435944, -0.00293951, 0.000505313, 0.011321, -0.00730551, -0.00701741, -0.00393573, 0.00568603, -0.00441037, 4.47696e-06, -0.00113014, -0.00816183, -0.0023439, -0.00880084, 0.0133834, -0.00747356, 0.000238353, -0.00703211, 0.0175569, -0.0147934, -0.00592084, -0.00104689, 0.00193023, -0.00584238, 0.0235858, -0.0226491, 0.0283542, -0.0117708, -0.00770908, 0.000631497, 0.0256777, -0.00827822, 0.00763801, -0.0114972, -0.0117758, 0.00620899, 0.000476455, 0.0098074, -0.0123425, 0.0061932, -0.00922346, -0.00122258, 0.00685395, -0.00947542, 0.00379197, -0.0111819, -0.0258622, -0.0163714, 0.000379023, 0.0148602, 0.00814239, -0.0123476, -0.00409123, -0.00435955, 0.00702159, 0.00772662, 0.00221592, -0.0133679, 0.000629221, 0.0107924, -0.00467678, 0.00886887, -0.00579663, -0.0129847, -0.000285545, -0.00291079, 0.0034805, -0.0125291, 0.00492484, 0.00943752, -0.00268436, 0.00437742, 0.000593377, 0.00267865, -0.0100421, 0.0205674, 0.0167491, -0.00362788, -0.00726209, 0.00809637, -0.00875982, 0.00461384, -0.00449473, -0.00354654, 0.00722712, 0.00852314, 0.0117073, -0.00692884, 0.0213975, 0.0334928, -0.00239757, -0.00799365, -0.012033, -0.0052173, -0.00640274, 0.00655552, 0.0154968, 0.00677434, 0.00371682, -0.0158272, -0.00131305, 0.00224898, -0.00613299, -8.1768e-05, 0.0014908, 0.000753673, 0.00937905, 0.00630095, -0.00683983, 0.00765506, -0.00235396, 0.00897751, -0.00446957, -0.0127854, 0.0117484, 0.00110111, -0.00572371, 0.0156879, -0.000479286, -0.00458104, -0.0274776, -0.00515632, 0.0135409, 0.00197834, -0.00629265, 0.0222731, 0.00483006, -0.0119772, -0.0102637, 0.00127444, 0.017825, -0.0117528, -0.00774486, 0.0111866, 0.00823291, 0.00232714, -0.00540678, -0.00386004, 0.018626, 0.00193769, -0.0170498, -0.000494642, -0.0109847, -0.00948904, 0.00920589, -0.0252527, -0.00757891, -0.00567621, 0.00145231, -0.00021339, 0.00920078, -0.00668544, 0.00988638, -0.0150407, -0.00280452, 0.000589543, -0.0157093, 0.00871256, -0.0127608, -7.94979e-05, -0.00127517, -0.000522599, 0.00459375, -0.00687795, -0.0100312, -0.00610632, 0.00741354, -0.00377502, -0.0256207, -0.0147448, 0.00208709, 0.00990315, -0.00960484, 0.00267547, 0.015975, -0.00591776, 0.0106281, 0.00938553, 0.00226961, 0.0104078, 0.0177395, -0.00325351, -0.00736686, 0.00421624, 0.0205931, 0.00180534, -0.00373851, -0.013874, -0.00314101, -0.0108416, -0.0129715, 0.00947178, -0.00383016, 0.00472404, 0.00660181, 0.00603737, -0.00712546, -0.0054734, -0.0112765, 0.00385673, 0.00522128, -0.00527631, -0.0199117, 0.0121419, 0.00194943, -0.00353265, -0.00268178, -0.00252151, -0.026003, 0.0126732, -0.00698022, 0.00279553, 0.00270233, -0.00144387, -0.00488562, -0.00379344, 0.000924234, 0.0144168, 0.00706693, 0.00138135, 0.00383296, -0.00195424, 0.00350079, -0.0072443, -0.00354873, 0.00484002, -0.00346443, 0.0134539, 0.00246892, 0.00691835, 5.77599e-05, 0.00442391, -0.00922533, 0.00247923, 0.0007008, -0.000262607, 0.0149447, -0.000940468, 0.0101989, 0.00267996, 0.0073331, 0.00194967, -0.00904243, -0.0109351, 0.0123939, -0.0254097, -0.0108378, 0.00221882, 0.0120624, -0.00982211, 0.00388979, -0.000570549, -0.00444891, -0.00844491, 0.000615206, 0.0084769, 0.0071135, -0.00787109, -0.00605049, -0.0221345, -0.0185275, 0.0052713, -0.0098071, -0.00508354, -0.00951543, -0.00416532, -0.00584327, -0.00834383, -0.0047703, -0.00443108, -0.000861755, 0.0145326, -0.00290669, -0.00473628, -0.00263845, -0.00754708, 0.00639671, 0.0131621, 0.0037658, -0.000861914, 0.0260894, 0.00282825, -0.0129467, -0.00709605, 0.00850824, 0.0211541, 0.0119274, -0.0010508, -0.00475818, 0.0166539, -0.0169803, -0.00975323, 0.00286802, -0.028938, -0.0090465, 0.000431254, 0.00142143, -0.00746636, 0.0175214, 0.00742451, -0.0107145, 0.00090967, -0.000118416, 0.0121911, -0.00809848, 0.00507367, 0.007608, -0.00310818, 0.0107207, -0.0156409, -0.0120762, -0.00714422, 0.0155586, -0.00496367, 0.00187297, -0.0111635, -0.00445046, -0.0108473, -0.00326216, -0.000790771, -0.00360563, -0.00113692, -0.00913437, 0.008472, -0.000685937, 0.00918652, 0.00845642, 0.0125608, 0.00832128, -0.00650871, -0.00115916, 0.00985428, -0.00581525, 0.00946452, -0.00223969, -0.0047589, 0.00193619, -0.00586358, -0.00238768, -0.00529425, 0.0120506, -0.00265635, -0.00128535, 0.0183106, -0.000547503, 0.00833686, -0.00196416, -0.0062162, -0.00179572, 0.0100507, -0.00390093, -0.0109185, -0.00239214, 0.00747579, -0.0233939, -0.00112812, -0.00307341, 0.00311576, -0.00279058, -0.0186793, -0.00769162, -0.00813319, -0.0116624, 0.00984687, -0.0185186, -0.00706547, 0.0106825, -0.00491416, -0.0116783, 0.0138808, 0.00916554, 0.00192349, 0.00675454, 0.015904, 0.00657191, -0.00516508, 0.00161628, -0.0185241, -0.00322557, 0.00638492, 0.000740254, 0.000103586, -0.0111372, 0.00200773, 0.00400462, 0.00758657, 0.00256225, -0.0123283, 0.012767, 0.00870453, -0.0125194, 0.019174, -0.00953392, 0.0106227, -0.00589873, 0.013007, 0.00953267, -0.0012614, 0.0020385, 0.0143346, -0.0161178, 0.0152188, -0.000271384, 0.00921705, -0.00244183, -0.00811168, -0.00940522, -0.00683866, 0.00624366, 0.00518851, 0.0133499, -0.000127415, -0.0242106, 0.00088361, -0.000993983, -0.0166327, 0.000209465, 0.00733893, -0.00894768, 0.0189583, -0.00532601, 0.00266947, 0.00635124, 0.00521534, -0.00443554, -0.0027043, 0.0118315, -0.00645307, 0.00275187, -0.00597412, 0.00293763, -0.00685451, 0.00670444, -0.00758654, -0.00171357, 0.00016647, 0.00415241, -0.0171463, 0.00710353, 0.0108248, -0.00352332, -0.0126008, -0.00445409, -0.00425396, 0.00469941, -0.0095242, -0.010252, 0.00201316, 0.00696668, -0.00944373, -0.00740156, -0.00601917, -0.00249543, -0.00550116, 0.00136189, 0.0060288, -0.0205721, -0.00957341, 0.00298704, 0.00743144, 0.0127754, 0.01174, 0.00760148, 0.00298727, 0.017184, 0.0152179, -0.000617601, 0.000863676, 0.00430346, -0.00205344, -0.00678251, 0.0149416, -0.00120295, -0.00467559, -0.017191, 0.00290194, -0.0228182, -0.00942073, 0.00243156, -0.00354242, -0.000236125, -0.00744653, -0.00305193, 0.0204858, -0.026756, 0.00355898, 0.0118392, 0.00078338, 0.00860617, 0.00164093, 0.00157132, 0.00864563, 0.00688557, 0.0156244, 0.00575976, 0.00115453, 0.00774618, 0.00224118, -0.00849043, 0.00728544, 0.0143132, -0.0116179, -0.0100805, 0.00991093, 0.0189471, 0.0118983, 0.00566126, 0.000102677, 0.00826605, 0.0208867, -0.0160984, 0.0193512, 0.0140348, 0.0145061, -0.00481359, -0.00688219, 9.30049e-05, -0.0128774, 0.00113972, 0.00369326, 0.0195469, -0.00332844, -0.0204236, -0.00337809, -0.00974712, 0.0119946, -0.00892576, 0.00186458, -0.00506792, -0.00834328, -0.00223765, 0.0210094, -0.000558188, -0.00189972, -0.0145841, -0.00209296, -0.00907228, -0.000276518, -0.0100514, -0.00924995, -0.0136941, -0.0084146, -0.0117754, 0.00079418, -0.00756051, 0.00461642, 0.00282106, 0.00120929, 0.00528317, 0.00267308, -0.00166291, 0.00903335, -0.00127763, 0.00825044, -0.0195938, -0.0051208, 0.000704956, -0.00673837, 0.0119819, 0.00274125, 0.00826672, 0.00877626, -0.00393697, 0.0111674, 0.0125706, -0.000520402, -0.00250833, 0.00211367, 0.00692461, 0.00283706, 0.00263649, 0.0036395, -0.00192564, 0.0182662, -0.0103054, 0.00238505, -0.000141803, 0.00194606, 0.00480954, 0.0195272, -0.0169332, -0.00250419, 0.00578147, -0.0195273, 0.00257959, -0.0046732, 0.00320296, -0.0100432, -0.0113119, -0.0172939, -0.00877084, -0.0172505, -0.000990046, 0.0146654, -0.00256396, -0.00172531, -0.00239401, 0.022928, 0.0183907, 0.00237446, -0.00592667, 0.00394362, -0.00951821, -0.00175899, 0.00759519, -0.0100693, 0.0187191, 0.0047909, 0.0124127, -0.00514528, 0.0179953, 0.00777841, -0.00572023, -0.00252904, 0.00118107, -0.00301326, -0.00946566, -0.000242634, -0.00195974, 0.00183929, -0.0193146, -0.0163441, -0.015359, -0.0149805, 0.00279401, -0.0142756, -0.00858574, 0.00965335, -0.0150562, -0.0147738, -0.00377309, 0.0042173, 0.00680105, 0.0177316, 0.00751812, 0.00438632, -0.00483578, 0.00888007, -0.000566399, 0.00233661, 0.00843893, 0.0168605, -0.00590584, 0.00157813, -0.0154866, 0.00553573, -0.00125404, -0.00450803, -0.000505595, -0.00178169, 0.00236804, -0.0113701, 0.000440195, 0.00325204, 0.00271326, 0.0103224, 0.00496743, 0.00503979, -0.00276522, 0.00944774, 0.00534092, 0.00715707, 0.00927121, 0.00779704, 0.0120519, -0.013669, 0.0031565, 0.0100765, 0.0134123, 0.00508031, 0.0114199, 0.00878356, 0.00684186, -0.00415671, -0.00332741, 0.00559858, 0.00916494, 0.0207005, -0.00534523, 0.0163499, 0.00469594, -0.0263033, 0.0120134, 0.00915616, 0.0195757, 0.00508856, 0.0141581, 0.00397566, -0.00816497, -0.011804, -0.00798932, -0.012321, 0.00370647, 0.0118071, -0.0067009, -0.00284561, 0.00561738, -0.00832706, 0.0259991, -0.00494587, -0.00742179, 0.00608499, -0.00948318, -0.0103867, -0.0108092, 0.0053772, -0.00721695, 0.00554144, 0.0111569, 0.000407877, -0.0087944, 0.00410716, 0.00759595, -0.00112604, 0.00981577, 0.00361152, -0.0199203, -0.0129385, 0.00239944, -0.0137231, -0.00672274, -0.0153511, -0.0051009, -0.0044967, 0.0164954, -0.00856568, -0.0150982, -0.00120542, 0.0148974, -0.00130997, -0.00360708, 0.000741691, 0.00331897, 0.0119874, 0.010928, -0.00608949, -0.00592497, -0.00996312, -0.0148426, -0.000130615, 0.0176637, 0.00382988, -0.00666377, 0.0106181, -0.00597681, -0.00883804, -0.0174931, 0.00957817, -0.015444, 0.0110085, 0.0066106, 0.00297802, 0.0176334, 0.0134688, 0.00817577, 0.00269194, -0.00786617, -0.00989893, 0.00144086, -0.00226464, 0.0166988, 0.000457771, 0.004719, -0.000891262, 0.00121418, 0.00680879, 0.00839736, -0.00200556, 0.00800338, -0.00320243, -0.00414775, 0.00108123, 0.0051069, 0.0131504, -0.00114125, 0.00483941, -0.00942051, 0.00784224, -0.00762005, 0.00294447, 0.0175732, -0.0134622, 0.00384521, -0.0114729, -0.00138036, -0.0143263, -0.0180066, -0.00843455, 0.0213703, 0.00307682, 0.0144778, -0.0061029, 0.00940555, 0.0129604, -0.0100221, -0.00238803, -0.000562867, -0.00401814, -0.00804801, -0.00416178, 0.0158892, -0.0107325};
    std::vector<float> default_data = {-0.00146382,0.0013453,-0.0187138,0.0046065,-0.00214253,0.00163712,-0.00827944,0.00298595,0.0105547,0.000102154,0.0117457,-0.00546841,-0.0104944,0.00660682,-0.00625276,0.0148596,-0.00829081,-0.0255912,-0.00888707,-0.00539781,0.0101922,-0.00628956,-0.00482589,0.00339587,-0.00121306,0.0210886,-0.00371003,-0.00287389,-0.0230144,-0.0105935,-0.000615274,0.0145502,0.0135433,0.00925328,-0.00243275,0.0151561,0.00197497,0.0100886,0.00439499,0.00438945,0.00645743,-0.00128149,-0.0168599,0.0177643,-0.00613857,0.00469861,-0.00582398,0.00668493,-0.00103692,0.00149386,0.00624049,0.0153727,0.0117067,0.0107825,-0.0205006,0.0117196,-0.0145473,0.00136395,-0.0111552,-0.0171463,0.0112422,-0.0173985,-0.0147975,-0.0158694,0.0148247,-0.00727862,0.00754843,-0.001128,0.00984235,0.00326633,-0.0103745,-0.000764704,-0.0208402,0.00389231,0.00243215,0.00455092,0.00275194,0.0291628,0.00272422,-0.0320464,0.0186225,-0.0209501,0.0105544,0.00310367,-1.22802e-05,0.00404831,-0.0108115,0.0141863,-0.00400148,0.00926096,-0.00358203,0.00126072,0.00387892,-0.00569566,-0.00634654,0.00882249,-0.00677104,0.00204175,0.0135715,-0.02453,-0.00315325,-0.00379922,-0.00608541,0.0135717,0.000195746,0.0132359,-0.000912438,-0.00208138,-0.0161209,0.00281664,0.00785215,-0.00316253,0.00353801,-0.00271609,-0.0177443,-0.000590157,-0.0153723,-0.000539041,0.00386642,0.00129153,0.0154308,-0.00388037,0.000763218,-0.00644921,-0.0168173,-0.00900313,-0.00383684,-0.00587587,-0.00130386,-0.0237398,0.0101014,-0.00981124,0.00465062,0.0117503,-0.0154914,-0.00944311,-0.00524505,0.0106121,0.0152626,-0.0108715,-0.00395067,-0.018694,0.00886654,-0.000245811,0.00281015,-0.0043059,-0.00282224,0.0098182,-0.00226851,0.00680064,0.00137044,-0.011258,-0.00171483,-0.000395065,0.0151453,-0.00126506,0.000306945,-0.00368253,0.00983281,-0.00149726,-0.00584821,0.17,0.0141966,-0.00127846,-0.00704493,0.0024827,0.0181503,0.01151,0.0132445,-0.0084942,0.0101566,0.0116719,-0.000821409,-0.000890962,-0.000626403,-0.00568668,0.00450628,-0.00503789,0.00950093,-0.00502095,0.0182362,-0.0093429,0.00246156,-0.0033775,0.00884862,-0.00787764,-0.00248331,-0.0038751,-0.00916308,-0.0128537,0.00527725,-0.013623,-0.0124964,-0.00396654,-0.00314504,-0.0119639,-0.0140579,0.00954706,-0.00564152,0.00435944,-0.00293951,0.000505313,0.011321,-0.00730551,-0.00701741,-0.00393573,0.00568603,-0.00441037,4.47696e-06,-0.00113014,-0.00816183,-0.0023439,-0.00880084,0.0133834,-0.00747356,0.000238353,-0.00703211,0.0175569,-0.0147934,-0.00592084,-0.00104689,0.00193023,-0.00584238,0.0235858,-0.0226491,0.0283542,-0.0117708,-0.00770908,0.000631497,0.0256777,-0.00827822,0.00763801,-0.0114972,-0.0117758,0.00620899,0.000476455,0.0098074,-0.0123425,0.0061932,-0.00922346,-0.00122258,0.00685395,-0.00947542,0.00379197,-0.0111819,-0.0258622,-0.0163714,0.000379023,0.0148602,0.00814239,-0.0123476,-0.00409123,-0.00435955,0.00702159,0.00772662,0.00221592,-0.0133679,0.000629221,0.0107924,-0.00467678,0.00886887,-0.00579663,-0.0129847,-0.000285545,-0.00291079,0.0034805,-0.0125291,0.00492484,0.00943752,-0.00268436,0.00437742,0.000593377,0.00267865,-0.0100421,0.0205674,0.0167491,-0.00362788,-0.00726209,0.00809637,-0.00875982,0.00461384,-0.00449473,0.16,0.00722712,0.00852314,0.0117073,-0.00692884,0.0213975,0.0334928,-0.00239757,-0.00799365,-0.012033,-0.0052173,-0.00640274,0.00655552,0.0154968,0.00677434,0.00371682,-0.0158272,0.03,0.05,-0.00613299,-8.1768e-05,0.0014908,0.000753673,0.00937905,0.00630095,-0.00683983,0.00765506,-0.00235396,0.00897751,-0.00446957,-0.0127854,0.0117484,0.00110111,-0.00572371,0.0156879,-0.000479286,-0.00458104,-0.0274776,-0.00515632,0.0135409,0.00197834,-0.00629265,0.0222731,0.00483006,-0.0119772,-0.0102637,0.00127444,0.017825,-0.0117528,-0.00774486,0.0111866,0.00823291,0.00232714,-0.00540678,-0.00386004,0.018626,0.00193769,-0.0170498,-0.000494642,-0.0109847,-0.00948904,0.00920589,-0.0252527,-0.00757891,-0.00567621,0.00145231,-0.00021339,0.00920078,-0.00668544,0.00988638,-0.0150407,-0.00280452,0.000589543,-0.0157093,0.00871256,-0.0127608,-7.94979e-05,-0.00127517,-0.000522599,0.00459375,-0.00687795,-0.0100312,-0.00610632,0.00741354,-0.00377502,-0.0256207,-0.0147448,0.00208709,0.00990315,-0.00960484,0.00267547,0.015975,-0.00591776,0.0106281,0.00938553,0.00226961,0.0104078,0.0177395,-0.00325351,-0.00736686,0.00421624,0.0205931,0.00180534,-0.00373851,-0.013874,-0.00314101,-0.0108416,-0.0129715,0.00947178,-0.00383016,0.00472404,0.00660181,0.00603737,-0.00712546,-0.0054734,-0.0112765,0.00385673,0.00522128,-0.00527631,-0.0199117,0.0121419,0.00194943,-0.00353265,-0.00268178,-0.00252151,-0.026003,0.0126732,-0.00698022,0.00279553,0.00270233,-0.00144387,-0.00488562,-0.00379344,0.000924234,0.0144168,0.00706693,0.00138135,0.00383296,-0.00195424,0.00350079,-0.0072443,-0.00354873,0.00484002,-0.00346443,0.0134539,0.00246892,0.00691835,5.77599e-05,0.00442391,-0.00922533,0.00247923,0.0007008,-0.000262607,0.0149447,-0.000940468,0.0101989,0.00267996,0.0073331,0.00194967,-0.00904243,-0.0109351,0.0123939,-0.0254097,-0.0108378,0.00221882,0.0120624,-0.00982211,0.00388979,-0.000570549,-0.00444891,-0.00844491,0.000615206,0.0084769,0.0071135,-0.00787109,-0.00605049,-0.0221345,-0.0185275,0.0052713,-0.0098071,-0.00508354,-0.00951543,-0.00416532,-0.00584327,-0.00834383,-0.0047703,-0.00443108,-0.000861755,0.0145326,-0.00290669,-0.00473628,-0.00263845,-0.00754708,0.00639671,0.0131621,0.0037658,-0.000861914,0.0260894,0.00282825,-0.0129467,-0.00709605,0.00850824,0.0211541,0.0119274,-0.0010508,-0.00475818,0.0166539,-0.0169803,-0.00975323,0.00286802,-0.028938,-0.0090465,0.000431254,0.00142143,-0.00746636,0.0175214,0.00742451,-0.0107145,0.00090967,-0.000118416,0.0121911,-0.00809848,0.00507367,0.007608,-0.00310818,0.0107207,-0.0156409,-0.0120762,-0.00714422,0.0155586,-0.00496367,0.00187297,-0.0111635,-0.00445046,-0.0108473,-0.00326216,-0.000790771,-0.00360563,-0.00113692,-0.00913437,0.008472,-0.000685937,0.00918652,0.00845642,0.0125608,0.00832128,-0.00650871,-0.00115916,0.00985428,-0.00581525,0.00946452,-0.00223969,-0.0047589,0.00193619,-0.00586358,-0.00238768,-0.00529425,0.0120506,-0.00265635,-0.00128535,0.0183106,-0.000547503,0.00833686,-0.00196416,-0.0062162,-0.00179572,0.0100507,0.04,-0.0109185,-0.00239214,0.00747579,-0.0233939,-0.00112812,-0.00307341,0.00311576,-0.00279058,-0.0186793,-0.00769162,-0.00813319,-0.0116624,0.00984687,-0.0185186,-0.00706547,0.0106825,-0.00491416,-0.0116783,0.0138808,0.00916554,0.00192349,0.00675454,0.015904,0.00657191,-0.00516508,0.00161628,-0.0185241,-0.00322557,0.00638492,0.000740254,0.000103586,-0.0111372,0.00200773,0.00400462,0.00758657,0.00256225,-0.0123283,0.012767,0.00870453,-0.0125194,0.019174,-0.00953392,0.0106227,-0.00589873,0.013007,0.00953267,-0.0012614,0.0020385,0.0143346,-0.0161178,0.0152188,-0.000271384,0.00921705,-0.00244183,-0.00811168,-0.00940522,-0.00683866,0.00624366,0.00518851,0.0133499,-0.000127415,-0.0242106,0.00088361,-0.000993983,-0.0166327,0.000209465,0.00733893,-0.00894768,0.0189583,-0.00532601,0.00266947,0.00635124,0.00521534,-0.00443554,-0.0027043,0.0118315,-0.00645307,0.00275187,-0.00597412,0.00293763,-0.00685451,0.00670444,-0.00758654,-0.00171357,0.00016647,0.00415241,-0.0171463,0.00710353,0.0108248,-0.00352332,-0.0126008,-0.00445409,-0.00425396,0.00469941,-0.0095242,-0.010252,0.00201316,0.00696668,-0.00944373,-0.00740156,-0.00601917,-0.00249543,-0.00550116,0.00136189,0.0060288,-0.0205721,-0.00957341,0.00298704,0.00743144,0.0127754,0.01174,0.00760148,0.00298727,0.017184,0.0152179,-0.000617601,0.000863676,0.00430346,-0.00205344,-0.00678251,0.0149416,-0.00120295,-0.00467559,-0.017191,0.00290194,-0.0228182,-0.00942073,0.00243156,-0.00354242,-0.000236125,-0.00744653,-0.00305193,0.0204858,-0.026756,0.00355898,0.0118392,0.00078338,0.00860617,0.00164093,0.00157132,0.00864563,0.00688557,0.0156244,0.00575976,0.00115453,0.00774618,0.00224118,-0.00849043,0.00728544,0.0143132,-0.0116179,-0.0100805,0.00991093,0.0189471,0.0118983,0.00566126,0.000102677,0.00826605,0.0208867,-0.0160984,0.0193512,0.0140348,0.0145061,-0.00481359,-0.00688219,9.30049e-05,-0.0128774,0.00113972,0.00369326,0.0195469,-0.00332844,-0.0204236,-0.00337809,-0.00974712,0.0119946,-0.00892576,0.00186458,-0.00506792,-0.00834328,-0.00223765,0.0210094,-0.000558188,-0.00189972,-0.0145841,-0.00209296,-0.00907228,-0.000276518,-0.0100514,-0.00924995,-0.0136941,-0.0084146,-0.0117754,0.00079418,0.04,0.00461642,0.00282106,0.00120929,0.00528317,0.00267308,-0.00166291,0.25,-0.00127763,0.00825044,-0.0195938,-0.0051208,0.000704956,-0.00673837,0.0119819,0.00274125,0.00826672,0.00877626,-0.00393697,0.0111674,0.0125706,-0.000520402,-0.00250833,0.00211367,0.00692461,0.00283706,0.00263649,0.0036395,-0.00192564,0.0182662,-0.0103054,0.00238505,-0.000141803,0.00194606,0.00480954,0.0195272,-0.0169332,-0.00250419,0.00578147,-0.0195273,0.00257959,-0.0046732,0.00320296,-0.0100432,-0.0113119,-0.0172939,-0.00877084,-0.0172505,-0.000990046,0.0146654,-0.00256396,-0.00172531,-0.00239401,0.022928,0.0183907,0.00237446,-0.00592667,0.00394362,-0.00951821,-0.00175899,0.00759519,-0.0100693,0.0187191,0.0047909,0.0124127,-0.00514528,0.0179953,0.00777841,-0.00572023,-0.00252904,0.00118107,-0.00301326,-0.00946566,-0.000242634,-0.00195974,0.00183929,-0.0193146,-0.0163441,-0.015359,-0.0149805,0.00279401,-0.0142756,-0.00858574,0.00965335,-0.0150562,-0.0147738,-0.00377309,0.0042173,0.00680105,0.0177316,0.00751812,0.00438632,-0.00483578,0.00888007,-0.000566399,0.00233661,0.00843893,0.0168605,-0.00590584,0.00157813,-0.0154866,0.00553573,-0.00125404,-0.00450803,-0.000505595,-0.00178169,0.00236804,-0.0113701,0.000440195,0.00325204,0.00271326,0.0103224,0.00496743,0.00503979,0.26,0.00944774,0.00534092,0.22,0.00927121,0.00779704,0.0120519,-0.013669,0.0031565,0.0100765,0.0134123,0.00508031,0.0114199,0.00878356,0.00684186,-0.00415671,-0.00332741,0.00559858,0.00916494,0.0207005,-0.00534523,0.0163499,0.00469594,-0.0263033,0.0120134,0.00915616,0.0195757,0.00508856,0.0141581,0.00397566,-0.00816497,-0.011804,-0.00798932,-0.012321,0.00370647,0.0118071,-0.0067009,-0.00284561,0.00561738,-0.00832706,0.0259991,-0.00494587,-0.00742179,0.00608499,-0.00948318,-0.0103867,-0.0108092,0.0053772,-0.00721695,0.00554144,0.0111569,0.000407877,-0.0087944,0.00410716,0.00759595,-0.00112604,0.00981577,0.00361152,-0.0199203,-0.0129385,0.00239944,-0.0137231,-0.00672274,-0.0153511,-0.0051009,-0.0044967,0.0164954,-0.00856568,-0.0150982,-0.00120542,0.0148974,-0.00130997,-0.00360708,0.000741691,0.00331897,0.0119874,0.010928,-0.00608949,-0.00592497,-0.00996312,-0.0148426,-0.000130615,0.0176637,0.00382988,-0.00666377,0.0106181,-0.00597681,-0.00883804,-0.0174931,0.00957817,-0.015444,0.0110085,0.0066106,0.00297802,0.0176334,0.0134688,0.00817577,0.00269194,-0.00786617,-0.00989893,0.00144086,-0.00226464,0.0166988,0.000457771,0.004719,-0.000891262,0.00121418,0.00680879,0.00839736,-0.00200556,0.00800338,-0.00320243,-0.00414775,0.00108123,0.0051069,0.0131504,-0.00114125,0.00483941,-0.00942051,0.00784224,-0.00762005,0.00294447,0.0175732,-0.0134622,0.00384521,-0.0114729,-0.00138036,-0.0143263,-0.0180066,-0.00843455,0.0213703,0.00307682,0.0144778,-0.0061029,0.00940555,0.0129604,0.26,-0.00238803,-0.000562867,-0.00401814,-0.00804801,-0.00416178,0.0158892,-0.0107325};
    std::vector<float> expected_from_data = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
};

/*----------------------------------------------------------------------------------------------*/

using namespace testing;

MATCHER_P(NearWithPrecision, precision, "") {
    return abs(get<0>(arg) - get<1>(arg)) < precision;
}

/*----------------------------------------------------------------------------------------------*/

TEST_F(DespikeTest,MacanaDespike){
    //Setup Vectors
    VecDoub hValues(nSamples,0.);
    VecDoub hSampleFlags(nSamples,1.);

    for(int i = 0; i<nSamples; i++) hValues[i] = default_data[i];

    //Get fake Data
    //MakeData<VecDoub>(hValues, nSamples, distmean, diststddev);

    //Add spikes
    //AddSpikes(hValues, nSamples, nSpikes, SpikeSigma, diststddev);

    //Macana Despike
    Macanadespike(hValues, hSampleFlags, SpikeLim, nSamples, samplerate, despikeWindow, isLowpassed, isDespiked);

    std::vector<double> hsf(&hSampleFlags[0], &hSampleFlags[0] + nSamples);

    EXPECT_THAT(hsf, Pointwise(NearWithPrecision(err), expected_from_data));

}

/*----------------------------------------------------------------------------------------------*/

TEST_F(DespikeTest,EigenDespike){
    //Setup Vectors
    Eigen::VectorXf hSampleFlags(nSamples);
    hSampleFlags.setOnes();

    Eigen::VectorXf hValues = Eigen::VectorXf::Map(default_data.data(), nSamples);

    //Get fake Data
    //MakeData<Eigen::VectorXf>(hValues, nSamples, distmean, diststddev);

    //Add spikes
    //AddSpikes(hValues, nSamples, nSpikes, SpikeSigma, diststddev);

    //Eigen Despike
    Eigendespike(hValues, hSampleFlags, SpikeLim, nSamples, samplerate, despikeWindow, timeConstant, isLowpassed, isDespiked);

    std::vector<double> hsf(&hSampleFlags[0], &hSampleFlags[0] + nSamples);

    EXPECT_THAT(hsf, Pointwise(NearWithPrecision(err), expected_from_data));
}

} // namespace
