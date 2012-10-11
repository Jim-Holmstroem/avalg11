#include "list.h"
#include <iostream>

namespace tsp
{
    void list::push_back(double val)
    {
       _nodes[_size].val = val;
       _nodes[_size].next = &_nodes[0];
       _nodes[0].prev = &_nodes[_size];

       if (_size != 0) {
           _nodes[_size - 1].next = &_nodes[_size];
           _nodes[_size].prev = &_nodes[_size - 1];
       }
       ++_size;
    }

    list::node* list::end() const 
    {
        if (_size == 0) {
            throw std::logic_error("list is of size 0");
        }

        return _nodes[0].prev;
    }

    list::node* list::begin() const 
    {
        if (_size == 0) {
            throw std::logic_error("list is of size 0");
        }

        return &_nodes[0];
    }

    list::node* list::get(size_t i) const
    {
        if (i >= _size) {
            throw std::logic_error("index out of bounds");
        }

        node *n = &_nodes[0];
        for(size_t j = 0; j < i; ++j, n = n->next);

        return n;
    }

    void list::two_opt_swap(node *a, node *b) const
    {
        node *c = a->next;
        node *d = b->next;

        c->prev = d; // no worries, will be taken care of in loop
        a->next = b;
        d->prev = c;

        node *current = b;
        node *prev = a;
        while(current != d) {
            // the swap
            current->next = current->prev;
            current->prev = prev;

            //the update
            prev = current;
            current = current->next;
        }
    }

    std::ostream& operator<<(std::ostream &out, const list &l)
    {
        list::node *n = l.begin();
        out << "FROM FRONT: " << std::endl;
        out << "(";
        for(size_t i = 0; i < l.size()-1; ++i) {
            out << n->val << ", ";
            n = n->next;
        }
        out << n->val << ")";

        out << std::endl;
        
        out << "FROM BACK: " << std::endl;
        n = l.end();
        out << "(";
        for(size_t i = 0; i < l.size()-1; ++i) {
            out << n->val << ", ";
            n = n->prev;
        }
        out << n->val << ")";
        return out;
    }
}
