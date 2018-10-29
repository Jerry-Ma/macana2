#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <random>

#include <CCfits/CCfits>

#include "mpfit.h"
#include "nr3.h"
#include "gaussFit.h"

namespace {

class GaussFitTest : public ::testing::Test
{
protected:
    GaussFitTest(): e2(sd()) {}
    ~GaussFitTest() override {}
    void SetUp() override {
    }
    void TearDown() override {}

    std::random_device sd;  // seed generator
    std::mt19937 e2;  // rand generator
    std::uniform_real_distribution<double> rand{-0.5, 0.5};
};

using namespace testing;

MATCHER_P(NearWithPrecision, precision, "") {
    return abs(get<0>(arg) - get<1>(arg)) < precision;
}

TEST_F(GaussFitTest, GaussFitMPFit) {
    size_t r = 9;
    size_t c = 9;
    size_t n = r * c;
    double s = 1.;  // pixel size
    std::vector<double> image(n, 0.);
    std::vector<double> unity(n, 1.);
    std::vector<double> sigma(n, 0.025);
    std::vector<double> residual(n);
    std::vector<double> ipos(n);
    std::vector<double> jpos(n);
    // setup coords center at (15., 10.)
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < c; ++j) {
            ipos[i * c + j] = s * (i - (r - 1) / 2.);
            jpos[i * c + j] = s * (j - (r - 1) / 2.);
        }
    struct vars_struct vars;
    vars.az = &ipos[0];
    vars.el = &jpos[0];
    vars.y = &image[0];
    vars.sigma = &unity[0];  // to produce a residual of negative gaussian
    // run the gauss funct with fwhm of 3.1 * s and 2.9 * s, centered at 0.1 and -0.1
    std::vector<double> params = {0., 1., 3.1 * s, 2.9 * s, 0.1 * s, -0.1 * s, 0.};
    std::vector<double> fitvar(7, 1.0);
    std::vector<double> fiterr(7, 0.);
    myfunct_gauss(n, 7, &params[0], &residual[0], nullptr, &vars);
    for (size_t i = 0; i < n; ++i) {
        image[i] -= residual[i] + sigma[0] * rand(e2);
    }
    std::vector<double> expected_image = {
       0.1600367986846205, 0.2609964730371924, 0.3321365121390888, 0.3896689679512127, 0.4175505641794933, 0.3788994885196257, 0.3181568598497313, 0.246343566180291 , 0.15632590237085  , 0.2499548825132084, 0.3611314109967023, 0.4938393866628206, 0.5763441383440523, 0.6086540569375619, 0.5635481091235378, 0.464084292843971 , 0.3452500058530137, 0.2341834466319097, 0.3123776885653153, 0.4941601642719758, 0.6421494176111386, 0.7662336387572405, 0.8064829805793318, 0.7333780076008504, 0.6127812681837885, 0.4603664557701133, 0.2810433505329943, 0.3735983827760868, 0.562628825214571 , 0.7617922740401993, 0.8872989644473327, 0.9353183248824375, 0.8833870533652972, 0.7225433231784979, 0.5351386232654979, 0.336497864440098 , 0.401632633309658 , 0.6185100787337061, 0.8043080043834396, 0.9648297821839503, 1.007387076188674 , 0.9253268860490997, 0.77707448014021  , 0.5703227231221681, 0.3655026524722184, 0.3992756999289704, 0.5859101897905105, 0.7727367372036791, 0.9059481535357726, 0.9601045679302548, 0.8968682059768152, 0.7414112261210883, 0.5509704202147953, 0.3457572982225425, 0.3299761801333522, 0.5012811442181527, 0.6671882526127019, 0.7976748581815311, 0.8302524756171115, 0.7803449281177793, 0.6292009064929912, 0.4758841746606695, 0.3107232214244371, 0.2593623135607939, 0.3924169231784841, 0.5249372916995235, 0.6102940432694947, 0.6472339407818956, 0.6077702175779408, 0.5062208832010623, 0.3588347527636233, 0.2370782982548471, 0.1746531098275094, 0.2697513553980692, 0.3571145599297929, 0.4282360120581097, 0.4601885823142943, 0.4261066486915479, 0.3541567030552099, 0.2658136272145991, 0.1583531916565155
    };
    EXPECT_THAT(image, Pointwise(NearWithPrecision(sigma[0]), expected_image));
    // do mpfit
    int status;
    mp_par pars[7];
    memset(&pars[0], 0, sizeof(pars));
    pars[1].limited[0] = 1;
    pars[1].limits[0] = 0.;
    pars[2].limited[0] = 1;
    pars[2].limits[0] = 0.;
    pars[3].limited[0] = 1;
    pars[3].limits[0] = 0.;
    // fix pa
    pars[6].fixed = 1;
    fitvar[6] = 0.;

    mp_result result;
    memset(&result, 0, sizeof(result));
    result.xerror = &fiterr[0];
    vars.sigma = &sigma[0];  // real sigma
    status = mpfit(&myfunct_gauss, n, 7, &fitvar[0], pars, nullptr, &vars, &result);
    if (status <= 0){
        cerr<<"GausFit(): Bad fit."<<endl;
    }
    //set the outputs
    for(size_t i = 0; i < n; ++i){
        fiterr[i] = result.xerror[i];
    }
    EXPECT_THAT(params, Pointwise(NearWithPrecision(sigma[0] * 2.), fitvar));
    for (size_t i = 0; i < 7; ++i) {
        std::cerr << "par[" << i << "] = " << fitvar[i] << " +/- " << fiterr[i] << std::endl;
    }
}


}  // namespace

