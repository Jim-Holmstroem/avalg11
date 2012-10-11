#include "point.h"

namespace tsp
{
    std::ostream& operator<<(std::ostream& out, point p)
    {
        out << "(" << p.x() << ", " << p.y() << ")"; 
        return out;
    }
}
