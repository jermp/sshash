#pragma once

#include "../util.hpp"

#include <vector>
#include <deque>

namespace sshash {

struct node {
    node()
        : id(constants::invalid_uint64)
        , front(constants::invalid_uint64)
        , back(constants::invalid_uint64)
        , sign(true)
        , chain_id(constants::invalid_uint64)
        , left(constants::invalid_uint64)
        , right(constants::invalid_uint64) {}

    node(uint64_t i, uint64_t f, uint64_t b, bool s = true)
        : id(i)
        , front(f)
        , back(b)
        , sign(s)
        , chain_id(constants::invalid_uint64)
        , left(constants::invalid_uint64)
        , right(constants::invalid_uint64) {}

    uint64_t id, front, back;

    bool sign;  // '+' --> forward
                // '-' --> backward

    uint64_t chain_id;
    uint64_t left, right;

    void print() const {
        std::cout << id << ":[" << front << "," << back << "," << (sign ? '+' : '-') << "]"
                  << std::endl;
    }
};

typedef std::deque<node> walk_t;
typedef std::vector<walk_t> walks_t;

}  // namespace sshash