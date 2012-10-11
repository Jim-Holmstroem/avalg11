#include "gtest/gtest.h"
#include "swaplist.h"

#include <iostream>

TEST(swaplist_test, trivial)
{
    tsp::swaplist swpl; 
    
    for (int i=0;i<8;++i)
        swpl.push_back(i);

    int i=0; 
    for (std::list<int>::iterator it=swpl.begin();it!=swpl.end();++it,++i)
        ASSERT_EQ(i,*it);
}

TEST(swaplist_test, simple_swap)
{
    tsp::swaplist swpl;
    swpl.push_back(0);
    swpl.push_back(1);
    swpl.push_back(5);
    swpl.push_back(4);
    swpl.push_back(3);
    swpl.push_back(2);
    swpl.push_back(6);
    swpl.push_back(7);

    std::list<int>::iterator from = swpl.begin(); std::advance(from,1);
    std::list<int>::iterator to   = swpl.begin(); std::advance(  to,6);
    tsp::swaplist::swap(from,to);
    

   // for (std::list<int>::iterator it=swpl.begin();it!=swpl.end();++it) {
   //     std::cout << *it;
   // } 
    
    int i=0;
    for (std::list<int>::iterator it=swpl.begin();it!=swpl.end();++it,++i)
        ASSERT_EQ(i,*it);
}

TEST(swaplist_test, simple_asymmtric_swap)
{
    tsp::swaplist swpl;
    swpl.push_back(1);
    swpl.push_back(5);
    swpl.push_back(4);
    swpl.push_back(3);
    swpl.push_back(2);
    swpl.push_back(6);
    swpl.push_back(7);
    swpl.push_back(0);
    
    std::list<int>::iterator from = swpl.begin(); std::advance(from,0);
    std::list<int>::iterator to   = swpl.begin(); std::advance(  to,5);
    tsp::swaplist::swap(from,to);
    
    int i=0;
    for (std::list<int>::iterator it=swpl.begin();it!=swpl.end();++it,++i)
        ASSERT_EQ((i+1)%8,*it);
}

TEST(swaplist_test,little_devil)
{
    tsp::swaplist swpl;
    swpl.push_back(0);
    swpl.push_back(2);
    swpl.push_back(1);
    swpl.push_back(3);
        
    std::list<int>::iterator from = swpl.begin(); std::advance(from,0);
    std::list<int>::iterator to   = swpl.begin(); std::advance(  to,3);
    tsp::swaplist::swap(from,to);

    int i=0; 
    for (std::list<int>::iterator it=swpl.begin();it!=swpl.end();++it,++i)
        ASSERT_EQ(i,*it);
}
