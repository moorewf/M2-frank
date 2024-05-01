
#include "aring-numberfield.hpp"
#include <gtest/gtest.h>
#include "ring.hpp"

TEST(NumberField, create)
{
  std::vector<long> f { 1, 0, 1};
  M2::ARingNumberField F(f);

  // EXPECT_EQ(ringName(F), "NumberField");
  // EXPECT_EQ(R.cardinality(), static_cast<size_t>(-1));
  // EXPECT_EQ(R.characteristic(), 0);
}
