#ifndef INCLUDE_LIST_H
#define INCLUDE_LIST_H

#include <stdexcept>
#include <iostream>
#include <cstdlib>

namespace tsp
{
    /* A circular double linked list */
    class list
    {
        public:
            struct node 
            {
                int val;
                node *prev;
                node *next;

            };

            list(size_t size) : _max_size(size), _size(0)
            {
                _nodes = (node *) calloc(sizeof(node), _max_size);
            }

            virtual ~list()
            {
                free(_nodes);
            }

            void push_back(double val);
            node* get(size_t i) const;
            node* end() const;
            node* begin() const;

            //helpers
            static void connect(node* a,node* b) { a->next=b; b->prev=a ;};

            size_t size() const { return _size; }

            void two_opt_swap(node *a, node *b) const;
        private:
            list(const list &other);
            list& operator=(const list& other);
            size_t _max_size;
            size_t _size;
            node *_nodes;
    };

    std::ostream& operator<<(std::ostream &out, const list &l);
}

#endif // INCLUDE_LIST_H
