#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <random>

#include <CCfits/CCfits>

#include "astron_utilities.h"
#include "Map.h"
#include "gaussFit.h"

namespace {

class MapTest : public ::testing::Test
{
protected:
    MapTest(): rand_e2(rand_rd()) {}
    ~MapTest() override {}
    void SetUp() override {
        mWt.assign(mNr, mNc, 1.);
        mRcp.resize(mNr);
        mCcp.resize(mNc);
        // populate the coords with an index sequence
        // center to be zero
        for (int i = 0; i < mNr; ++i) {
            mRcp[i] = mPixsz * (i - (mNr + 1.) / 2.);
        }
        for (int i = 0; i < mNc; ++i) {
            mCcp[i] = mPixsz * (i - (mNc + 1.) / 2.);
        }
        mMap = std::make_unique<Map>(mName, mNr, mNc, mPixsz, mWt, mRcp, mCcp);
    }
    void TearDown() override {}
    std::string mName = "test map";
    int mNr = 30;   // pix
    int mNc = 20;   // pix
    double mPixsz = 1.5 / 3600. * TWO_PI / 360.;  // 1.5 arcsec in radians
    MatDoub mWt;
    VecDoub mRcp;   // center 0, radians
    VecDoub mCcp;    // center 0, radians
    std::unique_ptr<Map> mMap;

    std::random_device rand_rd;
    std::mt19937 rand_e2;
    //std::knuth_b e2(rd());
    //std::default_random_engine e2(rd()) ;
    std::uniform_real_distribution<double> rand {-0.5, 0.5};

    static void dumpImage(const MatDoub& image, const std::string& fname)
    {
        long naxis    =   2;
        // naxes goes as x, y (j, i)
        std::vector<long> naxes(naxis);
        naxes[0] = image.ncols();
        naxes[1] = image.nrows();
        long nelems = naxes[0] * naxes[1];
        auto pFits = std::make_unique<CCfits::FITS>(fname, CCfits::Write);

        std::valarray<double> arr(nelems);
        CCfits::ExtHDU* imageExt = pFits->addImage(
                    "gaussian model", DOUBLE_IMG, naxes);
        for(int j=0;j<naxes[0];j++)
          for(int i=0;i<naxes[1];i++)
            arr[naxes[1]*i+j] = image[i][j];
        imageExt->write(1,nelems,arr);
    }
};

using namespace testing;
// Map Constructor
TEST_F(MapTest, MapFromData) {
    int n = mMap->getNPixels();
    int r = mMap->getNrows();
    int c = mMap->getNcols();
    double* im = mMap->image[0];
    double* wt = mMap->weight[0];
    // check image and weight are properly initialized
    EXPECT_THAT(std::vector<double>(im, im + n), Each(DoubleEq(0.0)));
    EXPECT_THAT(std::vector<double>(wt, wt + n), Each(DoubleEq(1.0)));
    // add an gaussian
    MatDoub rcp(r, c);  // az
    MatDoub ccp(r, c);  // el
    MatDoub res(r, c, 0.);  // store the returned offset, which is the negate of the guassian
    MatDoub unity(r, c, 1.); // make up a unity sigma matrix
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j) {
            rcp[i][j] = mMap->getRowCoordsPhys(i);
            ccp[i][j] = mMap->getColCoordsPhys(j);
        }
    struct vars_struct vars;
    vars.az = rcp[0];
    vars.el = ccp[0];
    vars.y = im;
    vars.sigma = unity[0]; // need this to make res being a negative gausian
    // dc offset, peak val, az fwhm (rad), el fwhm (rad), az off (rad), el off (rad), pa,
    double p[7] = {0., 1., 4.9 / 3600. * TWO_PI / 360., 5.1 / 3600. * TWO_PI / 360., 15.5 * mPixsz, 10.5 * mPixsz, 0.};
    // double p[7] = {0., 1., 4.9, 5.1, 1.5, 0.5, 0.};
    myfunct_gauss(n, 7, p, res[0], nullptr, &vars);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j) {
            mMap->image[i][j] -= res[i][j]; // + 0.025 * rand(rand_e2);
        }
    // dump the array
    dumpImage(mMap->image, "!test_map_gaussian_model.fits");
    // start fitting
    VecInt fixme(7, 0);
    VecDoub fixVals(7, 0.);
    VecDoub params(7);
    VecDoub params_err(7);
    mMap->weight.assign(n, c, 1. / (0.025 * 0.025));  // set to correct weight = 1/sigma^2
    // vars.sigma =mMap->weight[0];
    /*
    mpGaussFit(&params[0], &params_err[0], rcp[0], ccp[0], im, sigma[0], n,  fixme, fixVals, p);
    for (int i = 0; i < 7; ++i) {
        cout << "par[" << i << "] = " << params[i] << " +/- " << params_err[i + 7] << endl;
    }
    */
    // do the fit
    VecDoub pp(14, 0.);
    mMap->fitToGaussian(pp, fixme, fixVals, nullptr, 120);  // 120 arcsec enclose all 80 pix in radius
    // mMap->fitToGaussian(pp, fixme, fixVals, nullptr, 3600. * 100. * mMap->getPixelSize());
    for (size_t i = 0; i < 7; ++i) {
        cerr << "fitToGaussian: par[" << i << "] = " << pp[i] << " +/- " << pp[i + 7] << endl;
    }
    mMap->fitToGaussianMasked(pp, fixme, fixVals, nullptr, 120);  // 120 arcsec enclose all 80 pix in radius
    // mMap->fitToGaussian(pp, fixme, fixVals, nullptr, 3600. * 100. * mMap->getPixelSize());
    for (size_t i = 0; i < 7; ++i) {
        cerr << "fitToGaussianMasked: par[" << i << "] = " << pp[i] << " +/- " << pp[i + 7] << endl;
    }
}

// Tests that Foo does Xyz.
/*
TEST_F(AnalParamsTest, AnalParamsGetNFiles) {
  EXPECT_EQ(ap->getNFiles(), 1);
}
*/

}  // namespace

