#include "gtest/gtest.h"
#include "map.h"

TEST(map_test, test_distance)
{
    double x1 = 0, y1 = 0, x2 = 3, y2 = 4;
    tsp::map g(2);

    g.add_city(x1, y1);
    g.add_city(x2, y2);
    g.calculate_distances();

    ASSERT_EQ(5, g.distance(0, 1));
    ASSERT_EQ(5, g.distance(1, 0));
    ASSERT_EQ(0, g.distance(0, 0));
    ASSERT_EQ(0, g.distance(1, 1));
}

TEST(map_test, two_opt)
{
    tsp::map map(4);
    map.add_city(0, 0);
    map.add_city(4, 0);
    map.add_city(4, 3);
    map.add_city(0, 3);


    map.calculate_distances();

    std::list<int> path;
    path.push_back(0);
    path.push_back(2);
    path.push_back(3);
    path.push_back(1);

    ASSERT_EQ(18, map.path_length(path));

    map.two_opt(path);

    ASSERT_EQ(14, map.path_length(path));

}
