#pragma once

namespace sshash {

struct vertex {
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

    struct range {
        range() {}
        range(uint64_t b, uint64_t e, uint64_t p) : begin(b), end(e), position(p) {}
        uint64_t begin;
        uint64_t end;
        uint64_t position;
    };

    cover() : num_sequences(0) {}

    void compute(std::vector<vertex>& vertices,
                 uint64_t num_runs_abundances  // TODO: remove from here
    ) {
        essentials::timer_type timer;

        timer.start();

        /* (abundance, num_seqs_with_front_abundance=abundance) */
        std::unordered_map<uint64_t, range> abundance_map;
        std::vector<color> colors;
        std::vector<vertex> tmp_vertices;
        uint64_t num_runs = num_runs_abundances;
        std::cout << "initial number of runs = " << num_runs << std::endl;
        num_sequences = vertices.size();

        walk_t walk;
        walks_t walks_in_round;

        {
            /*
                optimize for the most frequent case:
                push vertices of the form v:[mfa,mfa] to the bottom, and remove them
                forming a single (usually, long) walk.
                Warning: here we are assuming the mfa is also the smallest abundance.
            */

            std::sort(vertices.begin(), vertices.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front > y.front;
                if (x.back != y.back) return x.back > y.back;
                return x.id > y.id;
            });

            uint64_t count = 0;
            uint64_t total = vertices.size();
            while (!vertices.empty()) {
                auto v = vertices.back();
                if (v.front == constants::most_frequent_abundance and
                    v.back == constants::most_frequent_abundance) {
                    walk.push_back(v);
                    vertices.pop_back();
                    count += 1;
                } else {
                    break;
                }
            }
            std::cout << "num vertices of the form v:[" << constants::most_frequent_abundance << ","
                      << constants::most_frequent_abundance << "] = " << count << "/" << total
                      << "(" << (count * 100.0) / total << "%)" << std::endl;

            /* create a new vertex for next round */
            tmp_vertices.emplace_back(walks_in_round.size(), walk.front().front, walk.back().back);
            walks_in_round.push_back(walk);

            walk.clear();
        }

        while (true) {
            std::cout << "round " << rounds.size() << std::endl;

            uint64_t num_vertices = vertices.size();
            std::cout << "  num_vertices " << num_vertices << std::endl;

            /* all unvisited */
            if (rounds.size() == 0) {
                /* remember: we removed some vertices but the id-space still spans
                   [0..num_sequences-1] */
                colors.resize(num_sequences);
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

            uint64_t prev_front = vertices.front().front;
            uint64_t prev_offset = 0;
            uint64_t offset = 0;
            for (auto const& x : vertices) {
                if (x.front != prev_front) {
                    abundance_map[prev_front] = {prev_offset, offset, 0};
                    prev_offset = offset;
                }
                offset += 1;
                prev_front = x.front;
            }
            abundance_map[prev_front] = {prev_offset, offset, 0};
            assert(offset == vertices.size());

            while (true) {
                walk.clear();

                /* take an unvisited vertex */
                uint64_t i = 0;
                for (; i != num_vertices; ++i) {
                    uint64_t id = vertices[i].id;
                    if (colors[id] == color::white) break;
                }

                if (i == num_vertices) break;  // all visited

                /* create a new walk */
                while (true) {
                    auto vertex = vertices[i];
                    uint64_t id = vertex.id;
                    assert(colors[id] != color::black);

                    colors[id] = color::gray;
                    uint64_t back = vertex.back;
                    if (walk.size() > 0) {
                        assert(walk.back().id != id);
                        colors[walk.back().id] = color::black;
                        colors[id] = color::black;
                    }
                    walk.push_back(vertex);

                    uint64_t offset = 0;
                    bool no_match_found = false;

                    auto it = abundance_map.find(back);

                    /* abundance back is not found, so no match is possible */
                    if (it == abundance_map.cend()) {
                        no_match_found = true;
                        break;
                    }
                    offset = (*it).second.begin + (*it).second.position;  // skip to position

                    /* search for a match */
                    while (true) {
                        if (offset == num_vertices) break;
                        auto candidate = vertices[offset];

                        /* skip */
                        if (candidate.id == id) {
                            offset += 1;
                            continue;
                        }

                        /* checked all candidate matches */
                        if (candidate.front != back) {
                            no_match_found = true;
                            break;
                        }

                        /* match found */
                        if (colors[candidate.id] != color::black) {
                            uint64_t& position = (*it).second.position;
                            assert((*it).second.begin + position <= offset);
                            position += 1;
                            break;
                        }

                        offset += 1;
                    }

                    assert(offset <= num_vertices);
                    if (no_match_found or offset == num_vertices) break;

                    /* valid match was found, then visit it next */
                    i = offset;
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
                num_runs -= num_mergings_in_round;
                std::cout << "  num_mergings = " << num_mergings_in_round << std::endl;
                std::cout << "  num_runs " << num_runs << std::endl;

                // std::cout << "created vertices in round " << rounds.size() << ":" << std::endl;
                // for (auto const& v : tmp_vertices) {
                //     std::cout << v.id << ":[" << v.front << "," << v.back << "]\n";
                // }
            }

            if (all_singletons) {
                std::cout << "STOP: all walks are singletons --> no new mergings were found"
                          << std::endl;
                break;
            }

            rounds.push_back(walks_in_round);

            vertices.swap(tmp_vertices);
            tmp_vertices.clear();
            walks_in_round.clear();
            abundance_map.clear();
        }

        timer.stop();

        std::cout << "cover computed in: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;

        // TODO: form final walks and check that all vertex id are present

        // uint64_t num_runs = num_runs_abundances;
        // std::cout << "computed " << walks.size() << " walks" << std::endl;
        // std::fill(colors.begin(), colors.end(), false);
        // for (auto const& walk : walks) {
        //     num_runs -= walk.size() - 1;
        //     uint64_t prev_back = walk.front().front;
        //     for (auto const& w : walk) {
        //         if (colors[w.id] == true) { std::cout << "ERROR: duplicate vertex." <<
        //         std::endl; } if (w.front != prev_back) { std::cout << "ERROR: path is broken." <<
        //         std::endl; } prev_back = w.back; colors[w.id] = true; std::cout << w.id << ":["
        //         << w.front << "," << w.back << "] ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "num_runs " << num_runs << std::endl;
    }

private:
    uint64_t num_sequences;
    std::vector<walks_t> rounds;
};

}  // namespace sshash