#ifndef TSP_PATH_H
#define TSP_PATH_H

#include<vector>
#include<algorithm>

namespace tsp {

    class path {
        public:
            
            path(size_t N) : _matrix(N,std::vector<int>(N)) {}

            void connect(int i,int j) {
               _matrix[i][j]=1;
               _matrix[j][i]=1;
            }
            void remove(int i,int j) {
               _matrix[i][j]=0;
               _matrix[j][i]=0;
                 
            }

            bool is_connected(int i,int j) {
                return _matrix[i][j];    
            }

            void swap(int i,int ip,int j,int jp) { //ii'->ij' and jj'->ji'
                remove(i,ip);
                remove(j,jp);
                connect(i,jp);
                connect(j,ip);
            }

            int reset_begin() {
                last = 0;
                return 0;
            }
           
            int begin() {
                return 0;    
            }

            int next() {
                int next;
                int choice1,choice2;
                
                //last is a little bit fugly

                for(;!_matrix[last][choice1];choice1++);
                choice2=choice1+1;
                for(;!_matrix[last][choice2];choice2++);
                
                if( choice1 == last ) {
                    last = at;
                    at = choice2;    
                } else {
                    last = at;
                    at = choice1;   
                }

                return next;
            }

            std::vector<int> get_path() {
                std::vector<int> path;
                path.reserve(_matrix.size());
                
                int at = reset_begin();
                
                do {
                    path.push_back(at); 
                    at = next(); 
                } while( at != begin() );

                return path;
            };

        private:
            std::vector< std::vector< int > > _matrix;
            int last;
    };
    
    
}

#endif
