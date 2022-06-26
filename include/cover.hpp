#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <deque>

#include "util.hpp"

namespace sshash {

namespace constants {
constexpr uint64_t most_frequent_weight = 1;  // usual value
}

struct node {
    // We assume there are less than 2^32 sequences and that
    // the largest weight fits into a 32-bit uint.
    node(uint32_t i, uint32_t f, uint32_t b, bool s) : id(i), front(f), back(b), sign(s) {}
    uint32_t id, front, back;
    bool sign;  // 'true'  --> '+'
                // 'false' --> '-'
};

struct cover {
    typedef std::deque<node> walk_t;
    typedef std::vector<walk_t> walks_t;

    cover(uint64_t num_sequences, uint64_t num_runs_weights)
        : m_num_sequences(num_sequences), m_num_runs_weights(num_runs_weights) {}

    void compute(std::vector<node>& nodes) {
        assert(nodes.size() == m_num_sequences);

        essentials::timer_type timer;
        timer.start();

        /* (weight, position of a candidate weight with front = weight) */
        std::unordered_map<uint32_t, uint32_t> weight_map;

        /* visited[node.id] = true if node has been visited; false otherwise */
        std::vector<bool> visited;

        /* set of unvisited nodes */
        std::unordered_set<uint32_t> unvisited_nodes;

        /* map from node id to offset+ into nodes */
        std::vector<uint32_t> id_to_offset;

        walk_t walk;
        unvisited_nodes.reserve(nodes.size());  // at most
        m_walks.reserve(nodes.size());          // at most

        std::cout << "initial number of runs = " << m_num_runs_weights << std::endl;

        /* add all the nodes with backward orientation */
        uint64_t num_nodes = nodes.size();
        for (uint64_t i = 0; i != num_nodes; ++i) {
            auto node = nodes[i];
            assert(node.sign == true);
            nodes.emplace_back(node.id, node.back, node.front, false);
        }

        num_nodes = nodes.size();
        std::cout << "  num_nodes " << num_nodes << std::endl;
        assert(num_nodes % 2 == 0);

        id_to_offset.resize(num_nodes / 2);
        visited.resize(num_nodes / 2);

        /* all nodes unvisited */
        std::fill(visited.begin(), visited.end(), false);

        for (auto const& node : nodes) unvisited_nodes.insert(node.id);
        assert(unvisited_nodes.size() == num_nodes / 2);

        std::sort(nodes.begin(), nodes.end(), [](auto const& x, auto const& y) {
            if (x.front != y.front) return x.front < y.front;
            if (x.back != y.back) return x.back < y.back;
            return x.id < y.id;
        });

        /* fill weight_map */
        {
            uint64_t prev_front = constants::invalid_uint64;
            uint64_t offset = 0;
            for (auto const& node : nodes) {
                if (node.front != prev_front) weight_map[node.front] = offset;
                offset += 1;
                prev_front = node.front;
            }
            assert(offset == nodes.size());
        }

        /* fill id_to_offset map */
        {
            uint64_t offset = 0;
            for (auto const& node : nodes) {
                if (node.sign) id_to_offset[node.id] = offset;
                offset += 1;
            }
        }

        uint64_t i = 0;  // position of an unvisited node in nodes

        /* some statistics */
        uint64_t total_nodes_visited_to_find_matches = 0;
        uint64_t total_nodes_visited_to_failure = 0;
        uint64_t total_nodes_visited_to_success = 0;
        uint64_t total_nodes = 0;

        while (!unvisited_nodes.empty()) {
            /* 1. take an unvisited node */
            {
                uint32_t unvisited_node_id = *(unvisited_nodes.begin());
                i = id_to_offset[unvisited_node_id];
            }

            /* 2. create a new walk */
            walk.clear();
            while (true) {
                total_nodes += 1;

                auto node = nodes[i];
                uint64_t id = node.id;
                assert(visited[id] == false);
                visited[id] = true;

                /* append node to walk */
                if (walk.empty()) {
                    walk.push_back(node);
                } else if (walk.back().back == node.front) {
                    walk.push_back(node);
                } else if (walk.back().back == node.back) {
                    change_orientation(node);
                    walk.push_back(node);
                } else if (walk.front().front == node.back) {
                    walk.push_front(node);
                } else if (walk.front().front == node.front) {
                    change_orientation(node);
                    walk.push_front(node);
                } else {
                    assert(false);
                }

                /* node has been visited, so erase it from unvisited_nodes */
                unvisited_nodes.erase(id);

                auto try_to_extend = [&](uint64_t back, uint64_t candidate_i) {
                    bool no_match_found = false;
                    uint64_t tmp_total_nodes_visited_to_failure = 0;
                    uint64_t tmp_total_nodes_visited_to_success = 0;

                    while (true) {
                        total_nodes_visited_to_find_matches += 1;
                        tmp_total_nodes_visited_to_failure += 1;
                        tmp_total_nodes_visited_to_success += 1;

                        if (candidate_i == num_nodes) {
                            total_nodes_visited_to_failure += tmp_total_nodes_visited_to_failure;
                            break;
                        }
                        auto candidate = nodes[candidate_i];

                        /* checked all candidate matches */
                        if (candidate.front != back) {
                            no_match_found = true;
                            total_nodes_visited_to_failure += tmp_total_nodes_visited_to_failure;
                            break;
                        }

                        /* skip */
                        if (candidate.id == id) {
                            candidate_i += 1;
                            continue;
                        }

                        /* match found */
                        if (visited[candidate.id] == false) {
                            total_nodes_visited_to_success += tmp_total_nodes_visited_to_success;
                            break;
                        }

                        candidate_i += 1;
                    }
                    assert(candidate_i <= num_nodes);

                    if (no_match_found or candidate_i == num_nodes) {
                        weight_map[back] = candidate_i;
                        return false;
                    }

                    /* valid match was found, then visit it next */
                    i = candidate_i;

                    /* update candidate position in weight_map to point to *next* position */
                    weight_map[back] = candidate_i + 1;

                    return true;
                };

                /* try to extend to the right */
                uint64_t candidate_i = weight_map[walk.back().back];
                bool found = try_to_extend(walk.back().back, candidate_i);

                /* if not possible, try to extend to the left */
                if (!found) {
                    candidate_i = weight_map[walk.front().front];
                    found = try_to_extend(walk.front().front, candidate_i);
                }

                if (!found) break;
            }
            assert(!walk.empty());
            m_walks.push_back(walk);
        }
        assert(unvisited_nodes.empty());

        if (total_nodes_visited_to_find_matches > total_nodes) {
            std::cout << "  avg. num. nodes visited to find a match = "
                      << static_cast<double>(total_nodes_visited_to_find_matches) / total_nodes
                      << std::endl;
            std::cout << "  avg. num. nodes visited to success = "
                      << static_cast<double>(total_nodes_visited_to_success) / total_nodes
                      << std::endl;
            std::cout << "  avg. num. nodes visited to failure = "
                      << static_cast<double>(total_nodes_visited_to_failure) / total_nodes
                      << std::endl;
        }

        std::cout << "  num_walks = " << m_walks.size() << std::endl;

        {
#ifndef NDEBUG
            std::fill(visited.begin(), visited.end(), false);
#endif
            uint64_t num_matches_in_round = 0;
            for (auto const& walk : m_walks) {
                num_matches_in_round += walk.size() - 1;
#ifndef NDEBUG
                uint64_t prev_back = walk.front().front;
                // std::cout << "=>";
                for (auto node : walk) {
                    if (visited[node.id]) std::cout << "ERROR: duplicate node." << std::endl;
                    if (node.front != prev_back) {
                        std::cout << "ERROR: path is broken." << std::endl;
                    }
                    prev_back = node.back;
                    visited[node.id] = true;
                    // std::cout << node.id << ":[" << node.front << "," << node.back << ","
                    //           << (node.sign ? '+' : '-') << "] ";
                }
                // std::cout << std::endl;
#endif
            }
            assert(m_num_runs_weights > num_matches_in_round);
            m_num_runs_weights -= num_matches_in_round;
            std::cout << "  num_matches = " << num_matches_in_round << std::endl;
            std::cout << "  num_runs " << m_num_runs_weights << std::endl;
        }

        timer.stop();
        std::cout << "cover computed in: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;
        std::cout << "final number of runs R = " << m_num_runs_weights << std::endl;

        assert(m_num_runs_weights >= 1);
    }

    void save(std::string const& filename) {
        std::ofstream out(filename.c_str());
        uint64_t num_sequences = 0;
        for (auto const& walk : m_walks) {
            for (auto node : walk) {
                // out << node.id << ":[" << node.front << "," << node.back << "] ";
                out << node.id;
                out << (node.sign ? " 1\n" : " 0\n");
            }
            num_sequences += walk.size();
        }
        if (num_sequences != m_num_sequences) {
            std::cerr << "Error: expected to write " << m_num_sequences << " but written "
                      << num_sequences << std::endl;
            throw std::runtime_error("wrong number of sequences written");
        }
        out.close();
    }

private:
    uint64_t m_num_sequences;
    uint64_t m_num_runs_weights;
    walks_t m_walks;

    void change_orientation(node& n) {
        std::swap(n.front, n.back);
        n.sign = !n.sign;
    }
};

}  // namespace sshash