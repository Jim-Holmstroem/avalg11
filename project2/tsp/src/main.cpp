#include <iostream>
#include "list.h"
#include "map.h"

int main(void)
{
    size_t n;
    std::cin >> n;
    tsp::map map(n);

    double x, y;
    for (size_t i = 0; i < n; ++i) {
        std::cin >> x >> y;
        map.add_city(x, y);
    }

    map.calculate_distances();
    
    tsp::list path(n);
    map.construct_path(path);
    map.optimize_path(path);

    tsp::list::node *node = path.begin();
    for(size_t i = 0; i < path.size(); ++i) {
       std::cout << node->val << std::endl;
       node = node->next;
    }

    return 0;
}
