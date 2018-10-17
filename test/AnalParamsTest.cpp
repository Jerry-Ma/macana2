#include "AnalParams.h"
#include "gtest/gtest.h"

namespace {

// The fixture for testing class Foo.
class AnalParamsTest : public ::testing::Test {
protected:
  AnalParamsTest() {
      ap = new AnalParams(apXmlFile);
  }
  ~AnalParamsTest() override {
      delete ap;
  }

  void SetUp() override {}
  void TearDown() override {}

  std::string apXmlFile = "test/data/test_ap.xml";
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

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
