#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../util.hpp"

namespace sshash {

struct node {
    // We assume there are less than 2^32 sequences and that
    // the largest abundance fits into a 32-bit uint.
    node(uint32_t i, uint32_t f, uint32_t b, bool s) : id(i), front(f), back(b), sign(s) {}
    uint32_t id, front, back;
    bool sign;  // '+' --> forward
                // '-' --> backward
};

struct cover {
    typedef std::vector<node> walk_t;
    typedef std::vector<walk_t> walks_t;

    cover(uint64_t num_sequences, uint64_t num_runs_abundances)
        : m_num_sequences(num_sequences), m_num_runs_abundances(num_runs_abundances) {}

    void compute(std::vector<node>& nodes) {
        assert(nodes.size() == m_num_sequences);

        essentials::timer_type timer;
        timer.start();

        /* (abundance, position of a candidate abundance with front = abundance) */
        std::unordered_map<uint32_t, uint32_t> abundance_map;

        /* visited[node.id] = true if node has been visited; false otherwise */
        std::vector<bool> visited;

        /* set of unvisited nodes */
        std::unordered_set<uint32_t> unvisited_nodes;

        /* map from node id to offset+ into nodes */
        std::vector<uint32_t> id_to_offset;

        /* nodes created for next round */
        std::vector<node> tmp_nodes;

        walk_t walk;
        walks_t walks_in_round;

        unvisited_nodes.reserve(nodes.size());     // at most
        walk.reserve(nodes.size());                // at most
        walks_in_round.reserve(2 * nodes.size());  // at most

        std::cout << "initial number of runs = " << m_num_runs_abundances << std::endl;

        {
            /*
                optimize for the most frequent case:
                push nodes of the form v:[mfa,mfa] to the bottom, and remove them
                forming a single (usually, long) walk.
                Warning: here we are assuming the mfa is also the *smallest* abundance.
            */

            essentials::timer_type timer;
            timer.start();

            std::sort(nodes.begin(), nodes.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front > y.front;
                if (x.back != y.back) return x.back > y.back;
                return x.id > y.id;
            });

            uint64_t num_special_nodes = 0;
            uint64_t num_nodes = nodes.size();
            while (!nodes.empty()) {
                auto v = nodes.back();
                if (v.front == constants::most_frequent_abundance and
                    v.back == constants::most_frequent_abundance) {
                    walk.push_back(v);
                    nodes.pop_back();
                    num_special_nodes += 1;
                } else {
                    break;
                }
            }
            std::cout << "num nodes of the form v:[" << constants::most_frequent_abundance << ","
                      << constants::most_frequent_abundance << "] = " << num_special_nodes << "/"
                      << num_nodes << "(" << (num_special_nodes * 100.0) / num_nodes << "%)"
                      << std::endl;

            if (num_special_nodes != 0) {
                /* create new nodes for next round */
                assert(!walk.empty());
                uint32_t id = walks_in_round.size();
                uint32_t front = walk.front().front;
                uint32_t back = walk.back().back;
                tmp_nodes.emplace_back(id, front, back, true);
                tmp_nodes.emplace_back(id, back, front, false);
                walks_in_round.push_back(walk);
                walk.clear();
            }

            /* add all the nodes with backward orientation */
            num_nodes = nodes.size();
            for (uint64_t i = 0; i != num_nodes; ++i) {
                auto node = nodes[i];
                assert(node.sign == true);
                nodes.emplace_back(node.id, node.back, node.front, false);
            }

            timer.stop();
            std::cout << "  time: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;
        }

        while (true) {
            std::cout << "round " << rounds.size() << std::endl;

            uint64_t num_nodes = nodes.size();
            std::cout << "  num_nodes " << num_nodes << std::endl;
            assert(num_nodes % 2 == 0);

            essentials::timer_type round_timer;  // total time of round
            round_timer.start();

            if (rounds.size() == 0) {
                /* remember: we removed some nodes but the id-space still spans
                   [0..m_num_sequences-1] */
                id_to_offset.resize(m_num_sequences);
                visited.resize(m_num_sequences);
            } else {
                id_to_offset.resize(num_nodes / 2);
                visited.resize(num_nodes / 2);
            }

            /* all nodes unvisited */
            std::fill(visited.begin(), visited.end(), false);

            for (auto const& node : nodes) unvisited_nodes.insert(node.id);
            assert(unvisited_nodes.size() == num_nodes / 2);

            std::sort(nodes.begin(), nodes.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front < y.front;
                if (x.back != y.back) return x.back < y.back;
                return x.id < y.id;
            });

            /* fill abundance_map */
            {
                uint64_t prev_front = constants::invalid;
                uint64_t offset = 0;
                for (auto const& node : nodes) {
                    if (node.front != prev_front) abundance_map[node.front] = offset;
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

            uint64_t i = 0;                        // position of an unvisited node in nodes
            bool no_more_matches_possible = true;  // to stop the algorithm

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
                    walk.push_back(node);

                    /* node has been visited, so erase it from unvisited_nodes */
                    if (unvisited_nodes.find(id) != unvisited_nodes.cend()) {
                        unvisited_nodes.erase(id);
                    }

                    auto search_match = [&](uint64_t back, uint64_t candidate_i) {
                        bool no_match_found = false;
                        uint64_t tmp_total_nodes_visited_to_failure = 0;
                        uint64_t tmp_total_nodes_visited_to_success = 0;

                        while (true) {
                            total_nodes_visited_to_find_matches += 1;
                            tmp_total_nodes_visited_to_failure += 1;
                            tmp_total_nodes_visited_to_success += 1;

                            if (candidate_i == num_nodes) {
                                total_nodes_visited_to_failure +=
                                    tmp_total_nodes_visited_to_failure;
                                break;
                            }
                            auto candidate = nodes[candidate_i];

                            /* skip */
                            if (candidate.id == id) {
                                candidate_i += 1;
                                continue;
                            }

                            /* checked all candidate matches */
                            if (candidate.front != back) {
                                no_match_found = true;
                                total_nodes_visited_to_failure +=
                                    tmp_total_nodes_visited_to_failure;
                                break;
                            }

                            /* match found */
                            if (visited[candidate.id] == false) {
                                total_nodes_visited_to_success +=
                                    tmp_total_nodes_visited_to_success;
                                break;
                            }

                            candidate_i += 1;
                        }
                        assert(candidate_i <= num_nodes);

                        /* update candidate position in abundance_map */
                        abundance_map[back] = candidate_i;

                        if (no_match_found or candidate_i == num_nodes) return false;

                        /* valid match was found, then visit it next */
                        i = candidate_i;
                        return true;
                    };

                    /* 3. search for a match */

                    /* first, try to match back abundance */
                    uint64_t back = node.back;
                    uint64_t candidate_i = abundance_map[back];
                    bool found = search_match(back, candidate_i);

                    /* if a match is not found and the walk is singleton,
                       then the only node in the walk could be matched on
                       front abundance */
                    if (!found and walk.size() == 1) {
                        back = node.front;
                        candidate_i = abundance_map[back];
                        found = search_match(back, candidate_i);
                        if (found) {
                            /* change orientation of the node */
                            walk[0] = {node.id, node.back, node.front, !node.sign};
                        }
                    }

                    if (!found) break;
                    no_more_matches_possible = false;
                }
                assert(!walk.empty());

                /* create new nodes for next round */
                uint32_t id = walks_in_round.size();
                uint32_t front = walk.front().front;
                uint32_t back = walk.back().back;
                tmp_nodes.emplace_back(id, front, back, true);
                tmp_nodes.emplace_back(id, back, front, false);
                walks_in_round.push_back(walk);
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

            std::cout << "  num_walks = " << walks_in_round.size() << std::endl;

            {
#ifndef NDEBUG
                std::fill(visited.begin(), visited.end(), false);
#endif
                uint64_t num_matches_in_round = 0;
                for (auto const& walk : walks_in_round) {
                    num_matches_in_round += walk.size() - 1;
#ifndef NDEBUG
                    uint64_t prev_back = walk.front().front;
                    std::cout << "=>";
                    for (auto const& w : walk) {
                        if (visited[w.id]) std::cout << "ERROR: duplicate node." << std::endl;
                        if (w.front != prev_back) {
                            std::cout << "ERROR: path is broken." << std::endl;
                        }
                        prev_back = w.back;
                        visited[w.id] = true;
                        std::cout << w.id << ":[" << w.front << "," << w.back << ","
                                  << (w.sign ? '+' : '-') << "] ";
                    }
                    std::cout << std::endl;
#endif
                }
                assert(m_num_runs_abundances > num_matches_in_round);
                m_num_runs_abundances -= num_matches_in_round;
                std::cout << "  num_matches = " << num_matches_in_round << std::endl;
                std::cout << "  num_runs " << m_num_runs_abundances << std::endl;

                round_timer.stop();
                std::cout << "  time: " << round_timer.elapsed() / 1000000 << " [sec]" << std::endl;

                // std::cout << "created nodes in round " << rounds.size() << ":" << std::endl;
                // for (auto const& v : tmp_nodes) {
                //     std::cout << v.id << ":[" << v.front << "," << v.back << ","
                //               << (v.sign ? '+' : '-') << "]\n";
                // }
            }

            rounds.push_back(walks_in_round);

            nodes.swap(tmp_nodes);
            tmp_nodes.clear();
            walks_in_round.clear();
            abundance_map.clear();

            if (no_more_matches_possible) {
                std::cout << "STOP: no more matches possible." << std::endl;
                break;
            }
        }

        timer.stop();
        std::cout << "cover computed in: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;
        std::cout << "final number of runs R = " << m_num_runs_abundances << std::endl;

        assert(m_num_runs_abundances >= 1);
    }

    void save(std::string const& filename) {
        std::ofstream out(filename.c_str());
        assert(rounds.size() > 0);
        int r = rounds.size() - 1;
        const auto& walks = rounds[r];
        uint64_t num_sequences = 0;
        for (uint64_t w = 0; w != walks.size(); ++w) {
            assert(walks[w].size() == 1);
            assert(walks[w].front().sign == true);
            num_sequences += visit(w, r, true, out);
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
    uint64_t m_num_runs_abundances;
    std::vector<walks_t> rounds;

    /* visit walk of index w in round of index r */
    uint64_t visit(int w, int r, bool sign, std::ofstream& out) const {
        if (r > 0) {
            assert(size_t(w) < rounds[r].size());
            auto const& walk = rounds[r][w];
            uint64_t num_sequences = 0;
            for (auto const& node : walk) {
                /*
                    + & + = +
                    - & - = +
                    + & - = -
                    - & + = -
                */
                bool new_sign = sign and node.sign;
                if (sign == false and node.sign == false) new_sign = true;
                num_sequences += visit(node.id, r - 1, new_sign, out);
            }
            return num_sequences;
        }
        /* print */
        assert(size_t(w) < rounds[0].size());
        auto const& walk = rounds[0][w];

        if (sign) {
            /* print forward */
            for (uint64_t i = 0; i != walk.size(); ++i) {
                auto const& node = walk[i];
                // out << node.id << ":[" << node.front << "," << node.back << "] ";
                out << node.id;
                out << (node.sign ? " 1\n" : " 0\n");
            }
        } else {
            /* print backward and change orientation */
            for (int64_t i = walk.size() - 1; i >= 0; --i) {
                auto const& node = walk[i];
                // out << node.id << ":[" << node.back << "," << node.front << "] ";
                out << node.id;
                out << (node.sign ? " 0\n" : " 1\n");
            }
        }

        return walk.size();
    }
};

}  // namespace sshash