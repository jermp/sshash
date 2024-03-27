#pragma once

#include "include/util.hpp"

#include <vector>
#include <deque>

namespace sshash {

namespace constants {
constexpr uint32_t invalid_uint32 = uint32_t(-1);
}

struct node {
    node()
        : id(constants::invalid_uint32)
        , front(constants::invalid_uint32)
        , back(constants::invalid_uint32)
        , sign(true)
        , chain_id(constants::invalid_uint32)
        , left(constants::invalid_uint32)
        , right(constants::invalid_uint32) {}

    node(uint32_t i, uint32_t f, uint32_t b, bool s = true)
        : id(i)
        , front(f)
        , back(b)
        , sign(s)
        , chain_id(constants::invalid_uint32)
        , left(constants::invalid_uint32)
        , right(constants::invalid_uint32) {}

    uint32_t id, front, back;

    bool sign;  // '+' --> forward
                // '-' --> backward

    uint32_t chain_id;
    uint32_t left, right;

    void print() const {
        std::cout << id << ":[" << front << "," << back << "," << (sign ? '+' : '-') << "]"
                  << std::endl;
    }
};

typedef std::deque<node> walk_t;
typedef std::vector<walk_t> walks_t;

}  // namespace sshash