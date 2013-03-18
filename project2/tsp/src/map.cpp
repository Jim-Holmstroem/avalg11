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

    // Running time: O(n^2)
    void map::nearest_neighbour(list &path)
    {
        size_t current_city = 0;
        std::vector<bool> visited(_size, false);
        size_t num_visited = 0;

        while (num_visited < _size) {
            num_visited++;
            visited[current_city] = true;
            path.push_back(current_city);

            double min = 1000000000; 
            size_t min_index = -1;
            for (size_t i = 0; i < _size; ++i) {
                double d = dist(current_city, i);
                if (d < min && d != 0 && !visited[i]) {
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
        three_opt(path);
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

    void map::three_opt(list &path)
    {

        long counter=0;
        bool better = true;
        long N = 330000000;
        while(better && counter < N) 
        {
            list::node *a=path.begin(); //TODO could be randomized
            list::node *b=a->next->next->next->next;
            list::node *c=a->next->next;
            better = false;
            for ( //make it wider
                int n=0;
                b->next!=a && !better && counter < N; 
                b=b->next,++n
                ) 
                {
                bool first = true;
                for ( //slide [a,b]
                    ;
                    (a != path.begin() || first) && !better && counter < N; //came around again
                    a=a->next,b = b->next //slide b alongside a
                    ) 
                    {
                       first = false;
                    for ( //look inside [a,b]
                        c=a->next->next;
                        ( c->next != b ) && !better && counter < N;
                        c=c->next,counter+=10
                        ) 
                        {
                            double edge0 = 
                                dist(a->val,a->next->val)+
                                dist(b->val,b->next->val)+
                                dist(c->val,c->next->val); 
                            double edge1 = 
                                dist(a->val,c->next->val)+
                                dist(c->val,b->next->val)+
                                dist(b->val,a->next->val);
                            double edge2 = 
                                dist(a->val,b->val)+
                                dist(c->val,b->next->val)+
                                dist(a->next->val,c->next->val);
                            
                            if ( edge1<edge0 && edge1<edge2 ) { 
                                //only legit in order a<c<b
                                list::node* an = a->next; //the original a->next
                                list::connect(a,c->next); //changes a->next
                                list::connect(c,b->next);
                                list::connect(b,an); //the old a->next
                                 
                                std::swap(b,c); 
                                better=true;
                                counter+=5;
                            } else if ( edge2 < edge0 ) {
                                list::node* an = a->next;
                                list::node* bn = b->next;
                                list::node* cn = c->next;
                                
                                list::node* current  = b->prev;
                                list::node* prev = b;
                                list::node* next;

                                while(current!=c) 
                                {
                                    next=current->prev; 
                                    list::connect(prev,current); 
                                    prev=current;
                                    current=next;
                                    counter+=4;
                                }

                                list::connect(a,b);
                                list::connect(cn,an);
                                list::connect(c,bn);
                                
                                better=true;
                                counter+=10;
                            }

                       }
                  } 
            }
        }
    }
    void map::two_opt(list &path)
    {
        bool improved_path = true; // only to enter while loop
        while(improved_path) {
            improved_path = false;
            list::node *a = path.begin();
            list::node *an = a->next;
            for (size_t i = 0; i < path.size(); ++i) {
                list::node *b = path.begin();
                list::node *bn = b->next;
                for (size_t j = 0; j < path.size(); ++j) {
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
