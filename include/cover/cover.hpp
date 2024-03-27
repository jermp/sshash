#pragma once

#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "include/util.hpp"
#include "node.hpp"
#include "even_frequency_weights.hpp"

namespace sshash {

struct cover {
    cover(uint64_t num_sequences, uint64_t num_runs_weights)
        : m_num_sequences(num_sequences), m_num_runs_weights(num_runs_weights) {
        assert(m_num_runs_weights >= m_num_sequences);
    }

    void swap(std::vector<node>& nodes) {
        m_nodes.swap(nodes);
        // m_nodes.reserve(m_nodes.size() * 2);
    }

    void compute() {
        assert(m_nodes.size() == m_num_sequences);

        std::cout << "initial number of runs = " << m_num_runs_weights << std::endl;
        std::cout << "num_nodes = " << m_nodes.size() << std::endl;

        essentials::timer_type timer;
        timer.start();
        pre_process();
        merge_even();
        greedy_cover();
        timer.stop();

        std::cout << "cover computed in: " << timer.elapsed() / 1000000 << " [sec] ("
                  << (timer.elapsed() * 1000) / m_num_sequences << " [ns/node])" << std::endl;
    }

    void save(std::string const& filename) {
        std::ofstream out(filename.c_str());
        uint64_t num_sequences = 0;
        for (auto& walk : m_walks) {
            uint32_t prev_back = walk.front().front;
            for (auto& u : walk) {
                if (u.chain_id != constants::invalid_uint32) {
                    num_sequences += save_chain(true, u, out, prev_back);
                } else if (u.left != constants::invalid_uint32 and
                           u.right != constants::invalid_uint32) {
                    num_sequences += save_tree(true, u, out, prev_back);
                } else {
                    num_sequences += save_leaf(u, out, prev_back);
                }
            }
        }
        if (num_sequences != m_num_sequences) {
            std::cout << "Error: expected to write " << m_num_sequences << " but written "
                      << num_sequences << std::endl;
            throw std::runtime_error("wrong number of sequences written");
        }
        out.close();

        m_num_runs_weights = m_num_runs_weights - m_num_sequences + m_walks.size();
        assert(m_num_runs_weights >= 1);
        std::cout << "final number of runs = " << m_num_runs_weights << std::endl;
    }

private:
    uint64_t m_num_sequences;
    uint64_t m_num_runs_weights;

    walks_t m_walks;   // final walks
    walks_t m_chains;  // chains of equal nodes

    std::vector<node> m_nodes;

    /* (w, set of ids of nodes where w appears as front or back) */
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> m_incidence;

    /* set of unvisited nodes */
    std::unordered_set<uint32_t> m_unvisited;

    void check_link_and_update(node const& u, uint32_t& prev_back) {
        if (u.front != prev_back) std::cout << "ERROR: path is broken." << std::endl;
        prev_back = u.back;
    }

    void change_orientation(node& u) {
        std::swap(u.front, u.back);
        u.sign = !u.sign;
    }

    void insert_node(node const& u, uint32_t offset) {
        m_unvisited.insert(offset);
        m_incidence[u.front].insert(offset);
        m_incidence[u.back].insert(offset);
    }

    void erase_node(node const& u, uint32_t offset) {
        m_unvisited.erase(offset);
        m_incidence[u.front].erase(offset);
        m_incidence[u.back].erase(offset);
    }

    uint64_t save_leaf(node const& u, std::ofstream& out, uint32_t& prev_back) {
        assert(u.left == constants::invalid_uint32 and u.right == constants::invalid_uint32);
        assert(u.chain_id == constants::invalid_uint32);
        check_link_and_update(u, prev_back);
        out << u.id;
        out << (u.sign ? " 1\n" : " 0\n");
        return 1;
    }

    uint64_t save_chain(bool parent_sign, node const& v, std::ofstream& out, uint32_t& prev_back) {
        assert(v.chain_id != constants::invalid_uint32);
        assert(v.chain_id < m_chains.size());
        auto& chain = m_chains[v.chain_id];
        bool new_sign = parent_sign == v.sign;
        if (new_sign) {
            for (auto const& u : chain) save_leaf(u, out, prev_back);
        } else {
            for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
                auto& u = *it;
                change_orientation(u);
                save_leaf(u, out, prev_back);
            }
        }
        return chain.size();
    }

    uint64_t save_tree(bool parent_sign, node& u, std::ofstream& out, uint32_t& prev_back) {
        if (u.left == constants::invalid_uint32 and u.right == constants::invalid_uint32) {  // leaf
            if (u.chain_id != constants::invalid_uint32) {
                return save_chain(parent_sign, u, out, prev_back);
            }
            if (parent_sign == false) change_orientation(u);
            return save_leaf(u, out, prev_back);
        }
        uint64_t left_subtree_size = 0;
        uint64_t right_subtree_size = 0;
        bool new_sign = parent_sign == u.sign;

        if (new_sign) {
            left_subtree_size = save_tree(new_sign, m_nodes[u.left], out, prev_back);
            right_subtree_size = save_tree(new_sign, m_nodes[u.right], out, prev_back);
        } else {
            right_subtree_size = save_tree(new_sign, m_nodes[u.right], out, prev_back);
            left_subtree_size = save_tree(new_sign, m_nodes[u.left], out, prev_back);
        }
        return left_subtree_size + right_subtree_size;
    }

    void pre_process() {
        std::vector<node> tmp;
        tmp.reserve(m_nodes.size());  // at most

        /*
            nodes (x,y) and (y,x) are considered as equal,
            hence, put all nodes in the form (x,y) where x <= y
        */
        for (auto& u : m_nodes) {
            if (u.front > u.back) change_orientation(u);
        }

        std::sort(m_nodes.begin(), m_nodes.end(), [](auto const& x, auto const& y) {
            if (x.front != y.front) return x.front < y.front;
            return x.back < y.back;
        });

        walk_t chain;
        uint32_t front = m_nodes.front().front;
        uint32_t back = m_nodes.front().back;

        node dummy;
        dummy.front = 0;
        dummy.back = 0;
        m_nodes.push_back(dummy);

        for (auto& u : m_nodes) {
            /* copy front and back */
            uint32_t u_front = u.front;
            uint32_t u_back = u.back;

            if (u.front != front or u.back != back) {
                assert(!chain.empty());
                if (chain.size() == 1) {
                    tmp.push_back(chain.front());
                } else if (front != back and chain.size() % 2 == 0) {
                    /* create two parent nodes */
                    node p1, p2;
                    p1 = chain.back();
                    if (chain.size() == 2) {
                        p2 = chain.front();
                    } else {
                        chain.pop_back();
                        p2.front = chain.front().front;
                        p2.back = chain.back().back;
                        p2.chain_id = m_chains.size();
                        m_chains.push_back(std::move(chain));
                    }
                    tmp.push_back(p1);
                    tmp.push_back(p2);
                } else {
                    /* create one parent node */
                    node p;
                    p.front = chain.front().front;
                    p.back = chain.back().back;
                    p.chain_id = m_chains.size();
                    tmp.push_back(p);
                    m_chains.push_back(std::move(chain));
                }
                chain.clear();
            }

            append_node_to_walk(u, chain);

            front = u_front;
            back = u_back;
        }

        std::cout << "num_chains = " << m_chains.size() << std::endl;

        m_nodes.swap(tmp);
        std::vector<node>().swap(tmp);

        /* fill m_unvisited and m_incidence */
        uint32_t offset_u = 0;
        for (auto const& u : m_nodes) {
            insert_node(u, offset_u);
            offset_u += 1;
        }

        /*
            merge nodes (w,w):
            first erase (w,w) and then merge it with another node from incidence[w]
        */
        offset_u = 0;
        for (auto& u : m_nodes) {
            if (u.front == u.back) {
                uint32_t w = u.front;
                auto const& incidence_w = m_incidence[w];
                if (incidence_w.size() == 1) {
                    offset_u += 1;
                    continue;
                }
                erase_node(u, offset_u);
                assert(incidence_w.size() >= 1);
                uint32_t offset_x = *incidence_w.begin();
                auto& x = m_nodes[offset_x];
                erase_node(x, offset_x);
                auto p = merge(x, u, w, offset_x, offset_u);
                uint32_t offset_p = m_nodes.size();
                m_nodes.push_back(p);
                insert_node(p, offset_p);
            }
            offset_u += 1;
        }
    }

    void merge_even() {
        even_frequency_weights efw;

        {
            std::unordered_map<uint32_t, uint32_t> freq;  // (w, frequency of w)
            uint32_t offset_u = 0;
            for (auto const& u : m_nodes) {
                if (m_unvisited.find(offset_u) != m_unvisited.cend()) {
                    if (auto it = freq.find(u.front); it == freq.cend()) {
                        freq[u.front] = 1;
                    } else {
                        (*it).second += 1;
                    }
                    if (auto it = freq.find(u.back); it == freq.cend()) {
                        freq[u.back] = 1;
                    } else {
                        (*it).second += 1;
                    }
                }
                offset_u += 1;
            }
            efw.build(freq);
        }

        while (efw.has_next()) {
            /* 1. take weight w of lowest even frequency */
            uint32_t w = efw.min();
            auto const& incidence_w = m_incidence[w];
            assert(!incidence_w.empty());

            /* there is one single node (w,w) in the connected component */
            if (incidence_w.size() == 1) continue;

            /* 2. take two nodes x and y from m_incidence[w] and merge them into a parent node p */
            assert(incidence_w.size() >= 2);
            auto it = incidence_w.begin();
            uint32_t offset_x = *it;
            ++it;
            uint32_t offset_y = *it;

            auto& x = m_nodes[offset_x];
            auto& y = m_nodes[offset_y];
            auto p = merge(x, y, w, offset_x, offset_y);
            erase_node(x, offset_x);
            erase_node(y, offset_y);

            uint32_t offset_p = m_nodes.size();
            m_nodes.push_back(p);

            /* 3. if parent node p is (ww,ww), then merge it with a node from m_incidence[ww] */
            if (p.front == p.back) {
                uint32_t ww = p.front;
                efw.decrease_freq(ww);
                auto const& incidence_ww = m_incidence[ww];
                if (!incidence_ww.empty()) {
                    uint32_t offset_xx = *(incidence_ww.begin());
                    auto& xx = m_nodes[offset_xx];
                    insert_node(p, offset_p);
                    uint32_t offset_yy = offset_p;
                    auto& yy = m_nodes[offset_p];
                    p = merge(xx, yy, ww, offset_xx, offset_yy);
                    erase_node(xx, offset_xx);
                    erase_node(yy, offset_yy);
                    offset_p = m_nodes.size();
                    m_nodes.push_back(p);
                }
            }

            insert_node(p, offset_p);
        }
        assert(efw.R[0].e == efw.T.size());

        /*
            4. at the end, for each connected component CC:
            - if CC has only even-freq. nodes, there is only a single node in the component,
              which cannot be merged, thus its weights will appear as endpoints;
            - otherwise, there are only nodes whose weights have odd frequency.
        */

        compute_lower_bound();
    }

    void greedy_cover() {
        walk_t walk;

        while (!m_unvisited.empty()) {
            /* 1. take an unvisited node's offset */
            uint32_t offset_u = *(m_unvisited.begin());

            /* 2. create a new walk */
            walk.clear();

            while (true) {
                auto u = m_nodes[offset_u];

                /* append the node to current walk and erase the node */
                append_node_to_walk(u, walk);
                erase_node(u, offset_u);

                auto try_to_extend = [&](uint32_t w) {
                    auto const& incidence_w = m_incidence[w];
                    if (incidence_w.empty()) return false;
                    offset_u = *(incidence_w.begin());
                    return true;
                };

                /* 3. search for a match */

                /* try to extend to the right */
                bool found = try_to_extend(walk.back().back);

                /* if not possible, try to extend to the left */
                if (!found) found = try_to_extend(walk.front().front);

                /* no extension is possible, hence a path has been created */
                if (!found) break;
            }
            assert(!walk.empty());
            m_walks.push_back(walk);
        }

        assert(m_unvisited.empty());

        std::cout << "num_walks = " << m_walks.size() << std::endl;
    }

    void append_node_to_walk(node& u, walk_t& walk) {
        /* clang-format off */
        if (walk.empty()) { walk.push_back(u); return; }
        if (walk.front().front == u.front or walk.back().back == u.back) change_orientation(u);
        if (walk.front().front == u.back) { walk.push_front(u); } else
        if (walk.back().back == u.front) { walk.push_back(u); }
        /* clang-format on */
    }

    /* merge nodes x and y on weight w and return a parent node p */
    node merge(node& x, node& y, uint32_t w, uint32_t offset_x, uint32_t offset_y) {
        if (x.front == w) change_orientation(x);
        if (y.back == w) change_orientation(y);
        node p;
        p.front = x.front;
        p.back = y.back;
        p.left = offset_x;
        p.right = offset_y;
        return p;
    }

    void compute_lower_bound() {
        struct info {
            uint32_t freq; /* freq of weight */

            /*
                If this flag is true, then w always appears in nodes of the
                form (w,w) and its freq is always even.
            */
            bool all_equal;
        };

        std::unordered_map<uint32_t, info> weights;

        uint32_t offset_u = 0;
        for (auto const& u : m_nodes) {
            if (m_unvisited.find(offset_u) != m_unvisited.cend()) {  // if not visited
                uint64_t front = u.front;
                uint64_t back = u.back;
                if (auto it = weights.find(front); it == weights.cend()) {
                    weights[front] = {1, true};
                } else {
                    (*it).second.freq += 1;
                }
                if (auto it = weights.find(back); it == weights.cend()) {
                    weights[back] = {1, true};
                } else {
                    (*it).second.freq += 1;
                }
                if (front != back) {
                    weights[back].all_equal = false;
                    weights[front].all_equal = false;
                }
            }
            offset_u += 1;
        }

        uint64_t num_endpoints = 0;
        for (auto const& p : weights) {
            /* Special case: weight only appears in nodes of the form (w,w), so count twice. */
            if (p.second.all_equal) {
                num_endpoints += 2;
                assert(m_incidence[p.first].size() == 1);
                continue;
            }
            /*
                If the excess of frequency is odd, then the weight will appear as
                end-point.
            */
            if (p.second.freq % 2 == 1) num_endpoints += 1;
        }
        assert(num_endpoints % 2 == 0);
        uint64_t num_walks = num_endpoints / 2;
        std::cout << "(estimated) num_walks = " << num_walks << std::endl;
    }
};

}  // namespace sshash