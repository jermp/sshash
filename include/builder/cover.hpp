#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../util.hpp"

namespace sshash {

struct vertex {
    // We assume there are less than 2^32 sequences and that
    // the largest abundance fits into a 32-bit uint.
    vertex(uint32_t i, uint32_t f, uint32_t b, bool s) : id(i), front(f), back(b), sign(s) {}
    uint32_t id, front, back;
    bool sign;  // '+' --> forward
                // '-' --> backward
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
        std::unordered_map<uint32_t, uint32_t> abundance_map;
        std::vector<color> colors;
        std::vector<vertex> tmp_vertices;
        std::unordered_set<uint32_t> unvisited_vertices;
        /* map from vertex id to (offset+,offset-)*/
        std::vector<std::pair<uint32_t, uint32_t>> id_to_offsets;
        walk_t walk;
        walks_t walks_in_round;

        unvisited_vertices.reserve(vertices.size());  // at most
        walk.reserve(vertices.size());                // at most
        walks_in_round.reserve(2 * vertices.size());  // at most

        std::cout << "initial number of runs = " << m_num_runs_abundances << std::endl;

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
                /* create new vertices for next round */
                assert(!walk.empty());
                uint32_t id = walks_in_round.size();
                uint32_t front = walk.front().front;
                uint32_t back = walk.back().back;
                tmp_vertices.emplace_back(id, front, back, true);
                tmp_vertices.emplace_back(id, back, front, false);
                walks_in_round.push_back(walk);
                walk.clear();
            }

            /* now  add all the vertices with backward orientation */
            num_vertices = vertices.size();
            for (uint64_t i = 0; i != num_vertices; ++i) {
                auto vertex = vertices[i];
                assert(vertex.sign == true);
                vertices.emplace_back(vertex.id, vertex.back, vertex.front, false);
            }

            timer.stop();
            std::cout << "  time: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;
        }

        while (true) {
            std::cout << "round " << rounds.size() << std::endl;

            uint64_t num_vertices = vertices.size();
            std::cout << "  num_vertices " << num_vertices << std::endl;
            assert(num_vertices % 2 == 0);

            essentials::timer_type round_timer;  // total time of round
            round_timer.start();

            /* all unvisited */
            if (rounds.size() == 0) {
                /* remember: we removed some vertices but the id-space still spans
                   [0..m_num_sequences-1] */
                id_to_offsets.resize(m_num_sequences);
                colors.resize(m_num_sequences);
                std::fill(colors.begin(), colors.end(), color::invalid);
                for (auto const& vertex : vertices) colors[vertex.id] = color::white;
            } else {
                id_to_offsets.resize(num_vertices / 2);
                colors.resize(num_vertices / 2);
                std::fill(colors.begin(), colors.end(), color::white);
            }

            for (auto const& vertex : vertices) unvisited_vertices.insert(vertex.id);
            std::cout << "  num_unvisited_vertices " << unvisited_vertices.size() << std::endl;
            assert(unvisited_vertices.size() == num_vertices / 2);

            std::sort(vertices.begin(), vertices.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front < y.front;
                if (x.back != y.back) return x.back < y.back;
                return x.id < y.id;
            });

            // if (num_vertices > 0)
            {
                uint64_t prev_front = constants::invalid;
                uint64_t offset = 0;
                for (auto const& vertex : vertices) {
                    if (vertex.front != prev_front) abundance_map[vertex.front] = offset;
                    offset += 1;
                    prev_front = vertex.front;
                }
                assert(offset == vertices.size());
            }

            /* fill id_to_offsets map */
            {
                uint64_t offset = 0;
                for (auto const& vertex : vertices) {
                    if (vertex.sign) {
                        id_to_offsets[vertex.id].first = offset;
                    } else {
                        id_to_offsets[vertex.id].second = offset;
                    }
                    offset += 1;
                }
            }

            uint64_t i = 0;  // position of an unvisited vertex in vertices

            uint64_t total_vertices_visited_to_find_matches = 0;
            uint64_t total_vertices_visited_to_failure = 0;
            uint64_t total_vertices_visited_to_success = 0;
            uint64_t total_vertices = 0;
            bool no_more_matches_possible = true;

            while (!unvisited_vertices.empty()) {
                /* 1. take an unvisited vertex */
                {
                    uint32_t unvisited_vertex_id = *(unvisited_vertices.begin());
                    i = id_to_offsets[unvisited_vertex_id].first;
                }

                /* 2. create a new walk */
                walk.clear();
                while (true) {
                    total_vertices += 1;

                    auto vertex = vertices[i];
                    uint64_t id = vertex.id;

                    // std::cout << "  visiting vertex " << id << " at offset " << i << std::endl;

                    assert(colors[id] != color::black);
                    colors[id] = color::gray;

                    /* vertex has been visited, so erase it from unvisited_vertices */
                    if (unvisited_vertices.find(id) != unvisited_vertices.cend()) {
                        unvisited_vertices.erase(id);
                    }

                    if (walk.size() > 0) {
                        assert(walk.back().id != id);
                        colors[walk.back().id] = color::black;
                        colors[id] = color::black;
                    }

                    walk.push_back(vertex);

                    auto search_match = [&](uint64_t back, uint64_t candidate_i) {
                        bool no_match_found = false;
                        uint64_t tmp_total_vertices_visited_to_failure = 0;
                        uint64_t tmp_total_vertices_visited_to_success = 0;

                        while (true) {
                            total_vertices_visited_to_find_matches += 1;
                            tmp_total_vertices_visited_to_failure += 1;
                            tmp_total_vertices_visited_to_success += 1;

                            if (candidate_i == num_vertices) {
                                total_vertices_visited_to_failure +=
                                    tmp_total_vertices_visited_to_failure;
                                break;
                            }
                            auto candidate = vertices[candidate_i];

                            /* skip */
                            if (candidate.id == id) {
                                candidate_i += 1;
                                continue;
                            }

                            /* checked all candidate matches */
                            if (candidate.front != back) {
                                no_match_found = true;
                                total_vertices_visited_to_failure +=
                                    tmp_total_vertices_visited_to_failure;
                                break;
                            }

                            /* match found */
                            if (colors[candidate.id] != color::black) {
                                total_vertices_visited_to_success +=
                                    tmp_total_vertices_visited_to_success;
                                break;
                            }

                            candidate_i += 1;
                        }

                        assert(candidate_i <= num_vertices);

                        constexpr uint64_t num_optimized_rounds_for_speed = 2;
                        if (rounds.size() <= num_optimized_rounds_for_speed) {
                            /* update candidate position in abundance_map */
                            abundance_map[back] = candidate_i;
                        }

                        if (no_match_found or candidate_i == num_vertices) return false;

                        /* valid match was found, then visit it next */
                        i = candidate_i;
                        return true;
                    };

                    /* 3. search for a match */
                    uint64_t back = vertex.back;
                    uint64_t candidate_i = abundance_map[back];
                    bool found = search_match(back, candidate_i);

                    if (!found and walk.size() == 1) {
                        back = vertex.front;
                        candidate_i = abundance_map[back];
                        found = search_match(back, candidate_i);
                        if (found) {
                            /* change orientation of the vertex */
                            walk[0] = {vertex.id, vertex.back, vertex.front, !vertex.sign};
                        }
                    }

                    if (!found) break;

                    no_more_matches_possible = false;
                }
                assert(!walk.empty());

                if (walk.size() == 1) {
                    assert(colors[walk.front().id] == color::gray);  // visited but not merged
                    continue;
                }

                /* invariant: all vertices belonging to a walk are all black */
                assert(std::all_of(walk.begin(), walk.end(),
                                   [&](vertex const& v) { return colors[v.id] == color::black; }));

                /* create new vertices for next round */
                uint32_t id = walks_in_round.size();
                uint32_t front = walk.front().front;
                uint32_t back = walk.back().back;
                tmp_vertices.emplace_back(id, front, back, true);
                tmp_vertices.emplace_back(id, back, front, false);
                walks_in_round.push_back(walk);
            }
            assert(unvisited_vertices.empty());

            if (total_vertices_visited_to_find_matches > total_vertices) {
                std::cout << "  avg. num. vertices visited to find a match = "
                          << static_cast<double>(total_vertices_visited_to_find_matches) /
                                 total_vertices
                          << std::endl;
                std::cout << "  avg. num. vertices visited to success = "
                          << static_cast<double>(total_vertices_visited_to_success) / total_vertices
                          << std::endl;
                std::cout << "  avg. num. vertices visited to failure = "
                          << static_cast<double>(total_vertices_visited_to_failure) / total_vertices
                          << std::endl;
            }

            /* add all gray vertices (singleton walks) */
            for (auto const& v : vertices) {
                if (colors[v.id] == color::gray) {
                    uint32_t id = walks_in_round.size();
                    uint32_t front = v.front;
                    uint32_t back = v.back;

                    if (v.sign) {
                        tmp_vertices.emplace_back(id, front, back, true);
                        tmp_vertices.emplace_back(id, back, front, false);
                    } else {
                        tmp_vertices.emplace_back(id, front, back, false);
                        tmp_vertices.emplace_back(id, back, front, true);
                    }

                    walk_t walk;
                    walk.push_back(v);
                    walks_in_round.push_back(walk);

                    colors[v.id] = color::black;
                }
            }

            std::cout << "  num_walks = " << walks_in_round.size() << std::endl;

            {
#ifndef NDEBUG
                std::fill(colors.begin(), colors.end(), color::white);
#endif
                uint64_t num_mergings_in_round = 0;
                for (auto const& walk : walks_in_round) {
                    num_mergings_in_round += walk.size() - 1;
#ifndef NDEBUG
                    uint64_t prev_back = walk.front().front;
                    std::cout << "=>";
                    for (auto const& w : walk) {
                        if (colors[w.id] == color::black) {
                            std::cout << "ERROR: duplicate vertex." << std::endl;
                        }
                        if (w.front != prev_back) {
                            std::cout << "ERROR: path is broken." << std::endl;
                        }
                        prev_back = w.back;
                        colors[w.id] = color::black;
                        std::cout << w.id << ":[" << w.front << "," << w.back << ","
                                  << (w.sign ? '+' : '-') << "] ";
                    }
                    std::cout << std::endl;
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
                //     std::cout << v.id << ":[" << v.front << "," << v.back << ","
                //               << (v.sign ? '+' : '-') << "]\n";
                // }
            }

            rounds.push_back(walks_in_round);

            vertices.swap(tmp_vertices);
            tmp_vertices.clear();
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
        for (uint64_t w = 0; w != walks.size(); ++w) num_sequences += visit(w, r, out);
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
    uint64_t visit(int w, int r, std::ofstream& out) const {
        if (r > 0) {
            assert(size_t(w) < rounds[r].size());
            auto const& walk = rounds[r][w];
            uint64_t num_sequences = 0;
            for (auto const& vertex : walk) num_sequences += visit(vertex.id, r - 1, out);
            return num_sequences;
        }
        /* print */
        assert(size_t(w) < rounds[0].size());
        auto const& walk = rounds[0][w];
        for (auto const& vertex : walk) {
            // out << vertex.id << ":[" << vertex.front << "," << vertex.back << "] ";
            out << vertex.id;
            out << (vertex.sign ? " +\n" : " -\n");
        }
        // out << '\n';
        return walk.size();
    }
};

}  // namespace sshash