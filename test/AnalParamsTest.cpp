#include <gtest/gtest.h>
// #include <gmock/gmock-matchers.h>

#include "AnalParams.h"

namespace {

class AnalParamsTest : public ::testing::Test
{
protected:
    AnalParamsTest(): ap(new AnalParams(apXmlFile, 1)) {}
    ~AnalParamsTest() override {delete ap;}
    void SetUp() override {}
    void TearDown() override {}
    std::string apXmlFile = "data/test_apb.xml";
    AnalParams* ap = nullptr;
};

TEST_F(AnalParamsTest, AnalParamsFromApXml) {
  EXPECT_EQ(ap->beammapping, 1);
}

TEST_F(AnalParamsTest, AnalParamsGetNFiles) {
  EXPECT_EQ(ap->getNFiles(), 1);
}

// more test for all the methods

}  // namespace
