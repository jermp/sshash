#pragma once

namespace sshash {

struct vertex {
    vertex(uint32_t s, uint32_t f, uint32_t b) : id(s), front(f), back(b) {}
    uint32_t id, front, back;
};

struct cover {
    typedef std::vector<vertex> walk_t;
    typedef std::vector<walk_t> walks_t;

    enum color_t {
        white = 0,  // unvisited
        gray = 1,   // visited by not merged
        black = 2   // visited and merged
    };

    cover() : num_sequences(0) {}

    void compute(std::vector<vertex>& vertices,
                 uint64_t num_runs_abundances  // TODO: remove from here
    ) {
        /* (abundance, num_seqs_with_front_abundance=abundance) */
        std::unordered_map<uint64_t, uint64_t> abundance_map;
        std::vector<color_t> colors;
        std::vector<vertex> tmp_vertices;
        uint64_t num_runs = num_runs_abundances;
        num_sequences = vertices.size();

        while (true) {
            std::cout << "round " << rounds.size() << std::endl;

            uint64_t num_vertices = vertices.size();
            std::cout << "num_vertices " << num_vertices << std::endl;
            tmp_vertices.clear();
            abundance_map.clear();

            /* all unvisited */
            colors.resize(num_vertices);
            std::fill(colors.begin(), colors.end(), color_t::white);

            std::sort(vertices.begin(), vertices.end(), [](auto const& x, auto const& y) {
                if (x.front != y.front) return x.front < y.front;
                if (x.back != y.back) return x.back < y.back;
                return x.id < y.id;
            });

            for (auto const& x : vertices) {
                auto it = abundance_map.find(x.front);
                if (it != abundance_map.cend()) {  // found
                    (*it).second += 1;
                } else {
                    abundance_map[x.front] = 1;
                }
            }

            {
                std::vector<std::pair<uint64_t, uint64_t>> offsets;  // (abundance,offset)
                offsets.reserve(abundance_map.size());
                for (auto const& p : abundance_map) { offsets.emplace_back(p.first, p.second); }
                assert(offsets.size() > 0);
                std::sort(offsets.begin(), offsets.end(),
                          [](auto const& x, auto const& y) { return x.first < y.first; });
                uint64_t offset = 0;
                for (auto const& p : offsets) {
                    uint64_t ab = p.first;
                    abundance_map[ab] = offset;
                    offset += p.second;
                }
            }

            walk_t walk;
            walks_t walks_in_round;

            while (true) {
                walk.clear();

                /* take an unvisited vertex */
                uint64_t i = 0;
                for (; i != num_vertices; ++i) {
                    uint64_t id = vertices[i].id;
                    if (colors[id] == color_t::white) break;
                }

                if (i == num_vertices) break;  // all visited

                /* create a new walk */
                while (true) {
                    uint64_t id = vertices[i].id;

                    // std::cout << "visiting vertex " << id << ":[" << vertices[i].front << ","
                    //           << vertices[i].back << "]" << std::endl;
                    // std::cout << "at offset " << i << std::endl;

                    if (colors[id] == color_t::black) break;

                    colors[id] = color_t::gray;
                    uint64_t front = vertices[i].front;
                    uint64_t back = vertices[i].back;
                    if (walk.size() > 0) {
                        assert(walk.back().id != id);
                        colors[walk.back().id] = color_t::black;
                        colors[id] = color_t::black;
                    }
                    walk.emplace_back(id, front, back);
                    // std::cout << "added " << id << " to walk-" << walks.size() << std::endl;

                    uint64_t offset = abundance_map[back];

                    /* search for a match */
                    bool no_match_found = false;
                    while (true) {
                        if (offset == num_vertices) break;
                        auto vertex = vertices[offset];
                        // std::cout << "  vertex " << vertex.id << ":[" << vertex.front << ","
                        //           << vertex.back << "]" << std::endl;

                        /* skip */
                        if (vertex.id == id) {
                            offset += 1;
                            continue;
                        }

                        /* checked all candidate matches */
                        if (vertex.front != back) {
                            // std::cout << "  MISfound match " << vertex.id << ":[" << vertex.front
                            //           << "," << vertex.back << "]" << std::endl;
                            // std::cout << "  at offset " << offset << std::endl;
                            // std::cout << "NO MATCHES FOUND" << std::endl;
                            no_match_found = true;
                            break;
                        }

                        /* match found */
                        if (colors[vertex.id] != color_t::black) {
                            // std::cout << "  found match " << vertex.id << ":[" << vertex.front
                            //           << "," << vertex.back << "]" << std::endl;
                            // std::cout << "  at offset " << offset << std::endl;
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
                    assert(colors[walk.front().id] == color_t::gray);  // visited but not merged
                    continue;
                }

                /* invariant: all vertices belonging to a walk are all black */
                assert(std::all_of(walk.begin(), walk.end(), [&](vertex const& v) {
                    return colors[v.id] == color_t::black;
                }));

                /* create a new vertex for next round */
                tmp_vertices.emplace_back(walks_in_round.size(), walk.front().front,
                                          walk.back().back);
                walks_in_round.push_back(walk);
            }

            /* add all gray vertices (singleton walks) */
            for (auto const& v : vertices) {
                if (colors[v.id] == color_t::gray) {
                    tmp_vertices.emplace_back(walks_in_round.size(), v.front, v.back);
                    walk_t walk;
                    walk.push_back(v);
                    walks_in_round.push_back(walk);
                }
            }

            std::cout << "num_walks_in_round " << rounds.size() << ": " << walks_in_round.size()
                      << std::endl;

            bool all_singletons = true;
            {
                std::fill(colors.begin(), colors.end(), color_t::white);
                for (auto const& walk : walks_in_round) {
                    if (walk.size() > 1) all_singletons = false;
                    num_runs -= walk.size() - 1;
                    uint64_t prev_back = walk.front().front;
                    std::cout << "=>";
                    for (auto const& w : walk) {
                        if (colors[w.id] == color_t::black) {
                            std::cout << "ERROR: duplicate vertex." << std::endl;
                        }
                        if (w.front != prev_back) {
                            std::cout << "ERROR: path is broken." << std::endl;
                        }
                        prev_back = w.back;
                        colors[w.id] = color_t::black;
                        std::cout << w.id << ":[" << w.front << "," << w.back << "] ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "num_runs " << num_runs << std::endl;

                std::cout << "created vertices in round " << rounds.size() << ":" << std::endl;
                for (auto const& v : tmp_vertices) {
                    std::cout << v.id << ":[" << v.front << "," << v.back << "]\n";
                }
            }

            if (all_singletons) {
                std::cout << "all walks are singletons: no new mergings were found" << std::endl;
                break;
            }

            rounds.push_back(walks_in_round);
            vertices.swap(tmp_vertices);
        }

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