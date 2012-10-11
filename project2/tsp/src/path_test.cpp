#include "gtest/gtest.h"
#include "path.h"

#include <vector>

TEST(path_test,trivial) {
    tsp::path p(8); 
    
    for (int i=0;i<8;++i)
        p.push_back(i);

    int i=0;
    int at = p.reset_begin();
    for ( ; i<8 ; ++i ) {
        ASSERT_EQ(i,at);
        at = p.next();
    }

    std::vector<int> ans;
    ans = p.get_path(); 
   
    i=0;
    for( std::vector<int>::iterator it=ans.begin(); 
         it!=ans.end();
         ++it,++i) {
        ASSERT_EQ(i,*it);
    }
}

