#include "map.h"

#include <set>
#include <algorithm>

#define MAX_ITER 2000

namespace tsp
{
    void map::add_city(double x, double y)
    {
        _cities.push_back(point(x, y));
    }

    void map::calculate_distances()
    {
        for (size_t i = 0; i < _size; ++i) {
            for (size_t j = 0; j < _size; ++j) {
                _distances[i][j] = i == j? 0 : _cities[i].distance(_cities[j]);
            }
        }
    }

    // Running time: O(n^2log(n))
    void map::nearest_neighbour(list &path)
    {
        size_t current_city = 0;
        std::set<int> visited;

        while (visited.size() < _size) { 
            visited.insert(current_city);
            path.push_back(current_city);

            double min = 1000000000; 
            size_t min_index = -1;
            for (size_t i = 0; i < _size; ++i) {
                double d = dist(current_city, i);
                if (d < min && d != 0 && 
                    visited.find(i) == visited.end()) {
                    min = d;
                    min_index = i;
                }
            }
            current_city = min_index;
        }
    }
   
    void map::construct_path(list &path)
    {
        nearest_neighbour(path);
    }

    void map::optimize_path(list &path)
    {
        two_opt(path);
    }

    double map::path_length(const list &path) const
    {
        double d = 0;
        list::node *n = path.begin();
        for (size_t i = 0; i < path.size(); ++i) {
            d += dist(n->val, n->next->val);
            n = n->next;
        }
        return d;
    }

    void map::two_opt(list &path)
    {
        bool improved_path = true; // only to enter while loop
        while(improved_path && _counter<MAX_ITER) {
            improved_path = false;
            list::node *a = path.begin();
            list::node *an = a->next;
            for (size_t i = 0; i < path.size() && _counter<MAX_ITER; ++i) {
                list::node *b = path.begin();
                list::node *bn = b->next;
                for (size_t j = 0; j < path.size() && _counter<MAX_ITER; ++j) {
                    if (a->val != b->val  && 
                        an->val != b->val && 
                        bn->val != a->val) {
                        double edge1 = dist(a->val, an->val) + 
                                       dist(b->val, bn->val);
                        double edge2 = dist(an->val, bn->val) + 
                                       dist(a->val, b->val);
                       
                       if (edge1 > edge2) {
                            path.two_opt_swap(a, b);
                            an = a->next;
                            bn = b->next;
                            improved_path = true;
                            ++_counter;
                        }
                    }
                    b = b->next;
                    bn = bn->next;
                }
                a = a->next;
                an = an->next;
            }
        }
    }

    std::ostream& operator<<(std::ostream& out, map g)
    {
        for (int i = 0; i < g.size(); ++i) {
            out << i << ": " << g.get_city(i);
            if (i != g.size() - 1)
                out << std::endl;
        }
        return out;
    }

}
