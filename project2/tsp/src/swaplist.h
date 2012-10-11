#ifndef INCLUDE_SWAPLIST_H
#define INCLUDE_SWAPLIST_H

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
namespace tsp {
    class swaplist
    {
    public:       
        swaplist() {
        };       
        void push_back(int i) {
            _list.push_back(i);    
        };
        static void swap(std::list<int>::iterator ik,std::list<int>::iterator jk)
        {
            ++ik;
            std::reverse(ik,jk); 
            --ik;
        };
        std::list<int>::iterator begin()
        {
            return _list.begin();    
        };
        std::list<int>::iterator end()
        {
            return _list.end();    
        };
    private:
        std::list<int> _list;
    };

}


/*                std::vector<int>::reverse_iterator rit;
                for( rit=_fwd.rbegin(); rit<_fwd.rend(); ++rit ) {
                    _bwd.push_back(*rit);
                }
*/
 
/*
            int get(size_t i) {
               return _fwd[i];  
            };

            int bget(size_t i) {
               return _bwd[i]; 
            };
*/

/*
            std::vector<int>::iterator bbegin()
            {
                return _bwd.begin();    
            }
            std::vector<int>::iterator bend()
            {
                return _bwd.end();
            }

            const std::vector<int>& get_vector() const{
                return _fwd;  
            };
*/
 
  //          std::vector<int> _fwd;
  //          std::vector<int> _bwd;
  
        /*std::vector<int>::iterator i;
                std::vector<int>::iterator j;
                if(ik<jk) {
                    i = ik+1;
                    j = jk;
                } else {
                    i = jk+1;
                    j = ik;
                }*/
/*
                std::vector<int>::iterator bi 
                = ( (_fwd.end()-1) - j)
                - _fwd.begin()
                + _bwd.begin();
                
                std::vector<int>::iterator bj 
                = ( (_fwd.end()-1) - i)
                - _fwd.begin()
                + _bwd.begin();
*/                

                //int tmp;
//                for ( ; i != j && bi != bj; ++i,++bi) {
                    //if(*i!=*bi) { //since we know all elements are unique
  //                   *i^=(*bi);
      //              *bi^=( *i);
    //                 *i^=(*bi);
                    //tmp = *bi;
                    //*bi = *i;
                    //*i = tmp;
                    
                    //}
        //        }
  

#endif //INCLUDE_SWAPLIST_H

