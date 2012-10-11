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
        //while(improved_path) //NOTE just running it once
       
     //   std::cout << path << std::endl;

        int counter=0;
        bool better = true;
        int N =10;
        while(better && counter < N) //++counter<100
        {
            list::node *a=path.begin(); //TODO could be randomized
            list::node *b=a->next->next->next->next;
            list::node *c=a->next->next;
            //std::cout << path << std::endl;   
            better = false;
            for ( //make it wider
                int n=0;
                n<2 && b->next!=a && !better && counter < N; //lower means more local 
                b=b->next,++n
                ) 
                {
              //  std::cout << "width:" << n << std::endl;
               // std::cout << "only twice" << std::endl;
                bool first = true;
                for ( //slide [a,b]
                    ;
                    (a != path.begin() || first) && !better && counter < N; //came around again
                    a=a->next,b = b->next //slide b alongside a
                    ) 
                    {
                //    std::cout << "slide to:" << a->val << std::endl;
                   // std::cout << ">slide" << std::endl;
                    //std::cout << n << std::endl;
                    for ( //look inside [a,b]
                        c=a->next->next;
                        ( c->next != b ) && !better && counter < N;
                        c=c->next,counter+=10
                        ) 
                        {
                   //     std::cout << "counter=" << counter << std::endl;
                 //          std::cout << "(a:" << a->val << ",a':" << a->next->val << "c:" << c->val << ",c':" << c->next->val << ",b:" << b->val << ",b':" << b->next->val << ")" << std::endl;
                            //std::cout << ">>pling" << std::endl;
                            double edge0 = 
                                dist(a->val,a->next->val)+
                                dist(b->val,b->next->val)+
                                dist(c->val,c->next->val); //original (1)(2)(3)=Id
                            double edge1 = 
                                dist(a->val,c->next->val)+
                                dist(c->val,b->next->val)+
                                dist(b->val,a->next->val); //1st perm (rotation perm (123))
                            double edge2 = 
                                dist(a->val,b->val)+
                                dist(c->val,b->next->val)+
                                dist(a->next->val,c->next->val); //2nd perm (what permutation is this?)
                            //check out why there is no 3rd perm

                            //TODO remember the b<->c switchup

                          //  std::cout << "(" << edge0 << "," << edge1 << "," << edge2 << ")" << std::endl;

                            if ( edge1<edge0 && edge1<edge2 ) { //three_opt swap, 1st (123)
                       //         std::cout << "(" << a->val << "," << b->val << "," << c->val << ")" << std::endl; 
                                
                                //only legit in order a<c<b
                                list::node* an = a->next; //the original a->next
                                list::connect(a,c->next); //changes a->next
                                list::connect(c,b->next);
                                list::connect(b,an); //the old a->next
                             //   std::cout << "great sucess for mother rusia" << std::endl;
                                 
                              // std::cout << "(" << a->val << "," << b->val << "," << c->val << ")" << std::endl; 
                              // std::cout << "switch" << std::endl;       
                             //   std::cout << path << std::endl;
            //the prevs will still be the same
                                //RESET the ordering (ABC)->(ACB) by renaming b<->c
                                std::swap(b,c); 
                              //  std::cout << "1st" << std::endl;

                     //           std::cout << "(" << a->val << "," << b->val << "," << c->val << ")" << std::endl; 
                     //           std::cout << path << std::endl;
                              //  std::cout << "swapdido" << std::endl;     
                             //   counter++;

                    //            std::cout << "1 done" << std::endl;
                                    better=true;
                                    counter+=5;
                            } else if ( edge2 < edge0 ) {
                                //2nd

                              //  std::cout << "2nd" << std::endl;
                                
                              //  std::cout << path << std::endl;
                              //  std::cout << "(" << a->val << "," << b->val << "," << c->val << ")" << std::endl; 
                                list::node* an = a->next;
                                list::node* bn = b->next;
                                list::node* cn = c->next;
                                
                                list::node* current  = b;
                                list::node* prev = b->prev;
                                list::node* prevprev;
                                
                                while(current!=c) //so it flips c' but not c
                                {
                                    //std::cout << current->val << std::endl;
                                    prevprev=prev->prev;
                                    list::connect(current,prev);
                                    prev=prevprev;
                                    current=prev;
                                    counter+=4;
                                }
                                
                                list::connect(a,b);
                                list::connect(cn,an);
                                list::connect(c,bn);


                               // std::cout << path << std::endl;
                      //          std::cout << "2 done" << std::endl;
                               // counter++;
                                better=true;
                                counter+=10;
                            }
                            //TODO make a reset suchas c resets and one slides it one step,instead of continuing at already checked thing,instead of continuing at already checked thingss

                        //    if(counter>100)
                        //        reset=true;
                     //   std::cout << "n=" << counter << std::endl;
                       }
                       first = false;
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
