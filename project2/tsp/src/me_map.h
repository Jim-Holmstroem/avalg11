#ifndef INCLUDE_MAP_H
#define INCLUDE_MAP_H

#include <vector>
#include <iostream>
#include "list.h"
#include "point.h"

namespace tsp 
{
    class map 
    { 
        public:
            map(size_t n) : _size(n), _distances(n, std::vector<double>(n))
            { _counter=0;}
            void add_city(double x, double y);
            void calculate_distances();
            void construct_path(list &path);
            void optimize_path(list &path);
            
            point get_city(int i) const { return _cities[i]; }
            int size() const { return _size; }
            double dist(int city1, int city2) const 
            {
                return _distances[city1][city2];
            }
            void two_opt(list &path);
            double path_length(const list &path) const;

        private:
            void nearest_neighbour(list &path);
            size_t _counter;//testing purpose only
            size_t _size;
            std::vector<point> _cities;
            std::vector<std::vector<double> > _distances;
    };

    std::ostream& operator<<(std::ostream& out, map g);
}

#endif // INCLUDE_MAP_H
