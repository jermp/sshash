#pragma once

namespace sshash {

struct vertex {
    vertex(uint32_t s, uint32_t f, uint32_t b) : seq_id(s), first_ab(f), last_ab(b) {}
    uint32_t seq_id, first_ab, last_ab;
};

struct cover {
    typedef std::vector<vertex> walk_t;

    void compute(std::vector<vertex>& vertices,
                 uint64_t num_runs_abundances  // TODO: remove from here
    ) {
        /* (abundance, num_seqs_with_front_abundance=abundance) */
        std::unordered_map<uint64_t, uint64_t> abundance_map;
        std::vector<bool> visited;

        while (true) {
            uint64_t num_vertices = vertices.size();

            std::sort(vertices.begin(), vertices.end(), [](auto const& x, auto const& y) {
                if (x.first_ab != y.first_ab) return x.first_ab < y.first_ab;
                if (x.last_ab != y.last_ab) return x.last_ab < y.last_ab;
                return x.seq_id < y.seq_id;
            });

            visited.resize(num_vertices);
            std::fill(visited.begin(), visited.end(), false);

            for (auto const& x : vertices) {
                auto it = abundance_map.find(x.first_ab);
                if (it != abundance_map.cend()) {  // found
                    (*it).second += 1;
                } else {
                    abundance_map[x.first_ab] = 1;
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

            while (true) {
                walk.clear();

                // take an unvisited vertex
                uint64_t i = 0;
                for (; i != num_vertices; ++i) {
                    uint64_t seq_id = vertices[i].seq_id;
                    if (visited[seq_id] == false) break;
                }

                if (i == num_vertices) break;  // all visited

                // std::cout << " i = " << i << std::endl;

                while (true) {
                    uint64_t seq_id = vertices[i].seq_id;
                    if (visited[seq_id] == true) break;

                    visited[seq_id] = true;
                    uint64_t front = vertices[i].first_ab;
                    uint64_t back = vertices[i].last_ab;
                    walk.emplace_back(seq_id, front, back);
                    // std::cout << "added " << seq_id << " to walk-" << walks.size() << std::endl;

                    uint64_t next_offset = abundance_map[back];

                    while (visited[seq_id] == true) {
                        next_offset += 1;
                        if (next_offset == num_vertices) break;
                        auto ab_endpoint = vertices[next_offset];
                        if (ab_endpoint.first_ab != back) break;
                        seq_id = ab_endpoint.seq_id;
                    }

                    if (next_offset == num_vertices) {
                        abundance_map[back] = next_offset - 1;
                        break;
                    }

                    if (vertices[next_offset].first_ab != back) {
                        abundance_map[back] = next_offset - 1;
                        i = next_offset - 1;
                    } else {
                        abundance_map[back] = next_offset;
                        i = next_offset;
                    }
                }

                if (!walk.empty()) {
                    // add walk to walks
                    walks.push_back(walk);
                    std::cout << "walks.size() = " << walks.size() << std::endl;
                }
            }

            break;
        }

        uint64_t num_runs = num_runs_abundances;
        std::cout << "computed " << walks.size() << " walks" << std::endl;
        std::fill(visited.begin(), visited.end(), false);
        for (auto const& walk : walks) {
            num_runs -= walk.size() - 1;
            uint64_t prev_back = walk.front().first_ab;
            for (auto const& w : walk) {
                if (visited[w.seq_id] == true) {
                    std::cout << "ERROR: duplicate vertex." << std::endl;
                }
                if (w.first_ab != prev_back) { std::cout << "ERROR: path is broken." << std::endl; }
                prev_back = w.last_ab;
                visited[w.seq_id] = true;
                std::cout << w.seq_id << ":[" << w.first_ab << "," << w.last_ab << "] ";
            }
            std::cout << std::endl;
        }
        std::cout << "num_runs " << num_runs << std::endl;
    }

private:
    std::vector<walk_t> walks;
};

}  // namespace sshash