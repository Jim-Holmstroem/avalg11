#include "map.h"

#include <set>
#include <algorithm>

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
    void map::nearest_neighbour(std::list<int> &path)
    {
        size_t current_city = 0;
        std::set<int> visited;

        while (visited.size() < _size) { 
            visited.insert(current_city);
            path.push_back(current_city);

            double min = 1000000000; 
            size_t min_index = -1;
            for (size_t i = 0; i < _size; ++i) {
                double d = distance(current_city, i);
                if (d < min && d != 0 && 
                    visited.find(i) == visited.end()) {
                    min = d;
                    min_index = i;
                }
            }
            current_city = min_index;
        }
    }
   
    void map::construct_path(std::list<int> &path)
    {
        nearest_neighbour(path);
    }

    void map::optimize_path(std::list<int> &path)
    {
        two_opt(path);
    }

    double map::path_length(const std::list<int> &path) const
    {
        double d = 0;
        std::list<int>::const_iterator now = path.begin();
        std::list<int>::const_iterator next= now;next++;

        for ( ; now != path.end() ; ++now,++next) {
            if(next==path.end())
                d += distance(*now,*path.begin()); //reconnect
            else
                d += distance(*now, *next);
            
        }
        return d;
    }

    void map::two_opt(std::list<int> &path)
    {
        if(path.size()<4)
            return;

        size_t counter=0;
        while(counter<(path.size()*path.size()+1)) 
        {

        std::list<int>::iterator i = path.begin();
        std::list<int>::iterator ip = i;
        ++ip; 
        for ( ;ip!=path.end();++i,++ip)
        {
            
            std::list<int>::iterator j = ip;
            std::advance(j,2);

            std::list<int>::iterator jp= j;
            jp--;

            for( ;j!=path.end();++j,++jp)
            {
                   // std::cout
                   // << std::distance(path.begin(),i) << "<"
                   // << std::distance(path.begin(),jp) << "<" 
                   // << std::distance(path.begin(),j) << "<" 
                   // << std::distance(path.begin(),path.end())
                   // << std::endl;
  
                 counter++;
                 //std::cerr << counter << std::endl;
                 if(distance(*i,*ip)+distance(*j,*jp) >
                     distance(*i,*jp)+distance(*j,*ip)
                    )
                  {
                    swaplist::swap(i,j);
//                    counter=0;
                //    break;
                  }
            }
             
        }
        }
        
        /*
        list::node *a = path.begin();
        list::node *ap = a->next;
        for (size_t i = 0; i < path.size(); ++i) {
            list::node *b = path.begin();
            list::node *bp = b->next;
            for (size_t j = 0; j < path.size(); ++j) {
                if (i == j || ap == b || bp == a) {
                    continue;
                }
                if (distance(a->val, ap->val) + distance(b->val, bp->val) >
                    distance(a->val, bp->val) + distance(b->val, ap->val)) {
                    std::cout << "path: " << path << std::endl;
                    std::cout << "a->val: " << a->val << " b->val: " << b->val << std::endl;
                    std::cout << "ap->val: " << ap->val << " bp->val: " << bp->val << std::endl;
                    path.swap(ap, bp);
                }
                b = b->next;
                bp = bp->next;
            }
            a = a->next;
            ap = ap->next;
        }
        */

      //  path = path.get_vector();

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
