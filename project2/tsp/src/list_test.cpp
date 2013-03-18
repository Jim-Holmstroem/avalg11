#include "gtest/gtest.h"
#include "list.h"

TEST(list_test, list_with_one_item)
{
    tsp::list list(1);
    list.push_back(2);

    ASSERT_EQ(2, list.begin()->val);
    ASSERT_EQ(2, list.get(0)->val);
    ASSERT_EQ(2, list.end()->val);
}

TEST(list_test, list_with_three_items)
{
    tsp::list list(3);
    list.push_back(1);
    list.push_back(2);
    list.push_back(3);

    ASSERT_EQ(1, list.begin()->val);
    ASSERT_EQ(1, list.get(0)->val);
    ASSERT_EQ(3, list.begin()->prev->val);
    ASSERT_EQ(2, list.get(1)->val);
    ASSERT_EQ(3, list.get(2)->val);
    ASSERT_EQ(3, list.end()->val);
    ASSERT_EQ(1, list.end()->next->val);
}
