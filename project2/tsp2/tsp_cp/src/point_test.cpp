#include "gtest/gtest.h"
#include "point.h"

TEST(point_test, create_point)
{
    double x = 4.4;
    double y = 5.5;

    tsp::point p(x, y);

    ASSERT_EQ(x, p.x());
    ASSERT_EQ(y, p.y());
}

TEST(point_test, test_distance)
{
    double x1 = 0, y1 = 0, x2 = 3, y2 = 4;
    tsp::point p1(x1, y1), p2(x2, y2);

    double distance = p1.distance(p2);

    ASSERT_EQ(distance, 5.0);
}
