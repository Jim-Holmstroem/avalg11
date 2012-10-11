#ifndef INCLUDE_MAP_H
#define INCLUDE_MAP_H

#include <vector>
#include <list>
#include <iostream>
#include "point.h"
#include "swaplist.h"

namespace tsp 
{
    class map 
    { 
        public:
            map(size_t n) : _size(n), _distances(n, std::vector<double>(n))
            { }
            void add_city(double x, double y);
            void calculate_distances();
            void construct_path(std::list<int> &path);
            void optimize_path(std::list<int> &path);
            
            point get_city(int i) const { return _cities[i]; }
            int size() const { return _size; }
            double distance(int city1, int city2) const 
            {
                return _distances[city1][city2];
            }
            void two_opt(std::list<int> &path);
            double path_length(const std::list<int> &path) const;

        private:
            void nearest_neighbour(std::list<int> &path);
            size_t _size;
            std::vector<point> _cities;
            std::vector<std::vector<double> > _distances;
    };

    std::ostream& operator<<(std::ostream& out, map g);
}

#endif // INCLUDE_MAP_H
