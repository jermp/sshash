#pragma once

#include <unordered_map>
#include <vector>

#include "../util.hpp"

namespace sshash {

struct vertex {
    // We assume there are less than 2^32 sequences and that
    // the largest abundance fits into a 32-bit uint.
    vertex(uint32_t s, uint32_t f, uint32_t b) : id(s), front(f), back(b) {}
    uint32_t id, front, back;
};

struct cover {
    typedef std::vector<vertex> walk_t;
    typedef std::vector<walk_t> walks_t;

    enum color {
        white = 0,  // unvisited
        gray = 1,   // visited by not merged
        black = 2,  // visited and merged
        invalid = -1
    };

    cover(uint64_t num_sequences, uint64_t num_runs_abundances)
        : m_num_sequences(num_sequences), m_num_runs_abundances(num_runs_abundances) {}

    void compute(std::vector<vertex>& vertices) {
        assert(vertices.size() == m_num_sequences);

        essentials::timer_type timer;

        timer.start();

        /* (abundance, position of a candidate abundance with front = abundance) */
        std::unordered_map<uint64_t, uint64_t> abundance_map;
        std::vector<color> colors;
        std::vector<vertex> tmp_vertices;
        std::cout << "initial number of runs = " << m_num_runs_abundances << std::endl;

        walk_t walk;
        walks_t walks_in_round;
        walk.reserve(vertices.size());            // at most
        walks_in_round.reserve(vertices.size());  // at most

        {
            /*
                optimize for the most frequent case:
                push vertices of the form v:[mfa,mfa] to the bottom, and remove them
                forming a single (usually, long) walk.
                Warning: here we are assuming the mfa is also the smallest abundance.
            */

            essentials::timer_type timer;
            timer.start();

            std::sort(vertices.begin(), vertices.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front > y.front;
                if (x.back != y.back) return x.back > y.back;
                return x.id > y.id;
            });

            uint64_t num_special_vertices = 0;
            uint64_t num_vertices = vertices.size();
            while (!vertices.empty()) {
                auto v = vertices.back();
                if (v.front == constants::most_frequent_abundance and
                    v.back == constants::most_frequent_abundance) {
                    walk.push_back(v);
                    vertices.pop_back();
                    num_special_vertices += 1;
                } else {
                    break;
                }
            }
            std::cout << "num vertices of the form v:[" << constants::most_frequent_abundance << ","
                      << constants::most_frequent_abundance << "] = " << num_special_vertices << "/"
                      << num_vertices << "(" << (num_special_vertices * 100.0) / num_vertices
                      << "%)" << std::endl;

            if (num_special_vertices != 0) {
                /* create a new vertex for next round */
                assert(!walk.empty());
                tmp_vertices.emplace_back(walks_in_round.size(), walk.front().front,
                                          walk.back().back);
                walks_in_round.push_back(walk);

                walk.clear();
            }

            timer.stop();
            std::cout << "  time: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;
        }

        while (true) {
            std::cout << "round " << rounds.size() << std::endl;

            uint64_t num_vertices = vertices.size();
            std::cout << "  num_vertices " << num_vertices << std::endl;

            essentials::timer_type round_timer;
            round_timer.start();

            /* all unvisited */
            if (rounds.size() == 0) {
                /* remember: we removed some vertices but the id-space still spans
                   [0..m_num_sequences-1] */
                colors.resize(m_num_sequences);
                std::fill(colors.begin(), colors.end(), color::invalid);
                for (auto const& v : vertices) colors[v.id] = color::white;
            } else {
                colors.resize(num_vertices);
                std::fill(colors.begin(), colors.end(), color::white);
            }

            std::sort(vertices.begin(), vertices.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front < y.front;
                if (x.back != y.back) return x.back < y.back;
                return x.id < y.id;
            });

            if (num_vertices > 0) {
                uint64_t prev_front = constants::invalid;
                uint64_t offset = 0;
                for (auto const& vertex : vertices) {
                    if (vertex.front != prev_front) abundance_map[vertex.front] = offset;
                    offset += 1;
                    prev_front = vertex.front;
                }
                assert(offset == vertices.size());
            }

            uint64_t i = 0;  // position of unvisited vertex in vertices

            while (true) {
                walk.clear();

                /* 1. take an unvisited vertex */

                // try to visit forward
                for (; i != num_vertices; ++i) {
                    uint64_t id = vertices[i].id;
                    if (colors[id] == color::white) break;
                }

                // if a vertex is not found, try to visit backward
                if (i == num_vertices) {
                    i = 0;
                    for (; i != num_vertices; ++i) {
                        uint64_t id = vertices[i].id;
                        if (colors[id] == color::white) break;
                    }
                }

                if (i == num_vertices) break;  // all visited

                /* 2. create a new walk */
                while (true) {
                    auto vertex = vertices[i];
                    uint64_t id = vertex.id;
                    assert(colors[id] != color::black);

                    colors[id] = color::gray;
                    if (walk.size() > 0) {
                        assert(walk.back().id != id);
                        colors[walk.back().id] = color::black;
                        colors[id] = color::black;
                    }
                    walk.push_back(vertex);

                    uint64_t back = vertex.back;
                    auto it = abundance_map.find(back);

                    /* abundance back is not found, so no match is possible */
                    if (it == abundance_map.cend()) break;

                    bool no_match_found = false;
                    uint64_t candidate_i = (*it).second;  // candidate position

                    /* 3. search for a match */
                    while (true) {
                        if (candidate_i == num_vertices) break;
                        auto candidate = vertices[candidate_i];

                        /* skip */
                        if (candidate.id == id) {
                            candidate_i += 1;
                            continue;
                        }

                        /* checked all candidate matches */
                        if (candidate.front != back) {
                            no_match_found = true;
                            break;
                        }

                        /* match found */
                        if (colors[candidate.id] != color::black) {
                            assert((*it).second <= candidate_i);
                            (*it).second += 1;
                            break;
                        }

                        candidate_i += 1;
                    }

                    assert(candidate_i <= num_vertices);
                    if (no_match_found or candidate_i == num_vertices) break;

                    /* valid match was found, then visit it next */
                    i = candidate_i;
                }

                assert(!walk.empty());

                if (walk.size() == 1) {
                    assert(colors[walk.front().id] == color::gray);  // visited but not merged
                    continue;
                }

                /* invariant: all vertices belonging to a walk are all black */
                assert(std::all_of(walk.begin(), walk.end(),
                                   [&](vertex const& v) { return colors[v.id] == color::black; }));

                /* create a new vertex for next round */
                tmp_vertices.emplace_back(walks_in_round.size(), walk.front().front,
                                          walk.back().back);
                walks_in_round.push_back(walk);
            }

            /* add all gray vertices (singleton walks) */
            for (auto const& v : vertices) {
                if (colors[v.id] == color::gray) {
                    tmp_vertices.emplace_back(walks_in_round.size(), v.front, v.back);
                    walk_t walk;
                    walk.push_back(v);
                    walks_in_round.push_back(walk);
                }
            }

            std::cout << "  num_walks = " << walks_in_round.size() << std::endl;

            bool all_singletons = true;
            {
#ifndef NDEBUG
                std::fill(colors.begin(), colors.end(), color::white);
#endif
                uint64_t num_mergings_in_round = 0;
                for (auto const& walk : walks_in_round) {
                    if (walk.size() > 1) all_singletons = false;
                    num_mergings_in_round += walk.size() - 1;
#ifndef NDEBUG
                    uint64_t prev_back = walk.front().front;
                    // std::cout << "=>";
                    for (auto const& w : walk) {
                        if (colors[w.id] == color::black) {
                            std::cout << "ERROR: duplicate vertex." << std::endl;
                        }
                        if (w.front != prev_back) {
                            std::cout << "ERROR: path is broken." << std::endl;
                        }
                        prev_back = w.back;
                        colors[w.id] = color::black;
                        // std::cout << w.id << ":[" << w.front << "," << w.back << "] ";
                    }
                    // std::cout << std::endl;
#endif
                }
                assert(m_num_runs_abundances > num_mergings_in_round);
                m_num_runs_abundances -= num_mergings_in_round;
                std::cout << "  num_mergings = " << num_mergings_in_round << std::endl;
                std::cout << "  num_runs " << m_num_runs_abundances << std::endl;

                round_timer.stop();
                std::cout << "  time: " << round_timer.elapsed() / 1000000 << " [sec]" << std::endl;
                // std::cout << "created vertices in round " << rounds.size() << ":" << std::endl;
                // for (auto const& v : tmp_vertices) {
                //     std::cout << v.id << ":[" << v.front << "," << v.back << "]\n";
                // }
            }

            rounds.push_back(walks_in_round);

            vertices.swap(tmp_vertices);
            tmp_vertices.clear();
            walks_in_round.clear();
            abundance_map.clear();

            if (all_singletons and rounds.size() > 1) {
                std::cout << "STOP: all walks are singletons --> no new mergings were found"
                          << std::endl;
                break;
            }
        }

        timer.stop();
        std::cout << "cover computed in: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;

        assert(m_num_runs_abundances >= 1);
    }

    void save(std::string const& filename) {
        std::ofstream out(filename.c_str());
        assert(rounds.size() > 0);
        int r = rounds.size() - 1;
        const auto& walks = rounds[r];
        for (auto const& walk : walks) {
            for (auto const& vertex : walk) visit(vertex.id, r, out);
        }
        out.close();
    }

private:
    uint64_t m_num_sequences;
    uint64_t m_num_runs_abundances;
    std::vector<walks_t> rounds;

    /* visit walk of index w in round of index r */
    void visit(int w, int r, std::ofstream& out) const {
        if (r > 0) {
            assert(size_t(w) < rounds[r].size());
            auto const& walk = rounds[r][w];
            for (auto const& vertex : walk) { visit(vertex.id, r - 1, out); }
        } else {  // print
            assert(size_t(w) < rounds[0].size());
            auto const& walk = rounds[0][w];
            for (auto const& vertex : walk) {
                // out << vertex.id << ":[" << vertex.front << "," << vertex.back << "] ";
                out << vertex.id << '\n';
            }
            // out << '\n';
        }
    }
};

}  // namespace sshash