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

// AnalParams Constructor
TEST_F(AnalParamsTest, AnalParamsFromApXml) {
  EXPECT_EQ(ap->beammapping, 0);
}

// Tests that Foo does Xyz.
TEST_F(AnalParamsTest, AnalParamsGetNFiles) {
  EXPECT_EQ(ap->getNFiles(), 1);
}

TEST(testcase, testset)
{
    EXPECT_EQ(1, 1);
}

}  // namespace
