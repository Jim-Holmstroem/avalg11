#ifndef INCLUDE_POINT_H
#define INCLUDE_POINT_H

#include <iostream>
#include <cmath>

namespace tsp
{
    class point
    {
        public:
            point() : _x(0), _y(0) {}
            point(double x, double y) : _x(x), _y(y) {}
            double x() const { return _x; }
            double y() const { return _y; }
            double distance(const point &p) const
            {
                double dx = _x - p.x();
                double dy = _y - p.y();
                return sqrt(pow(dx, 2) + pow(dy, 2));
            }
        private:
            double _x;
            double _y;
    };

    std::ostream& operator<<(std::ostream& out, point p);
}

#endif // INCLUDE_POINT_H
