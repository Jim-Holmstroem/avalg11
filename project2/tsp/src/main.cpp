#include <vector>
#include <iostream>
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
    
    std::list<int> path;
    map.construct_path(path);
    map.optimize_path(path);
    

    std::cerr << map.path_length(path) << std::endl;

    for(std::list<int>::iterator it = path.begin();it!=path.end();++it)
        std::cout << *it << std::endl;

    //tsp::list::node *node = path.begin();
    //for(size_t i = 0; i < path.size(); ++i) {
    //    std::cout << node->val << std::endl;
    //    node = node->next;
    //}

    return 0;
}
