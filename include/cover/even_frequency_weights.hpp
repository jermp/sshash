#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>  // for std::sort
#include <cassert>

namespace sshash {

struct even_frequency_weights {
    even_frequency_weights() {}

    struct wf {
        uint32_t w, f;
    };
    struct range {
        uint32_t b, e;
    };

    void build(std::unordered_map<uint32_t, uint32_t> const& freq) {
        uint32_t num_distinct_weights = freq.size();
        T.reserve(num_distinct_weights);
        P.reserve(num_distinct_weights);
        R[0] = {0, 0};

        for (auto const& p : freq) {
            if (p.second % 2 == 0) T.push_back({p.first, p.second});  // even
        }

        std::sort(T.begin(), T.end(), [](auto const& x, auto const& y) { return x.f < y.f; });

        uint32_t f = T.front().f;
        range r;
        r.b = 0;
        r.e = 0;
        for (uint32_t i = 0; i != T.size(); ++i) {
            P[T[i].w] = i;
            if (T[i].f == f) {
                r.e += 1;
            } else {
                R[f] = r;
                r.b = i;
                r.e = i + 1;  // one-past the end
            }
            f = T[i].f;
        }
        R[f] = r;  // last one

        // print();
    }

    void print() const {
        for (auto p : T) { std::cout << "(" << p.w << "," << p.f << ")"; }
        std::cout << std::endl;
        for (auto p : R) {
            std::cout << "R[" << p.first << "] = [" << p.second.b << "," << p.second.e << "]"
                      << std::endl;
        }
        for (auto p : P) { std::cout << "P[" << p.first << "] = " << p.second << std::endl; }
    }

    bool has_next() { return R[0].e < T.size(); }

    uint32_t min() {
        assert(has_next());
        uint32_t w = T[R[0].e].w;
        decrease_freq(w);
        return w;
    }

    void decrease_freq(uint32_t w) {
        uint32_t i = R[0].e;  // pos of weight of minimum frequency

        if (w != T[i].w) {
            if (P.find(w) == P.cend()) return;
            uint32_t j = P[w];
            assert(T[j].w == w);
            uint32_t f = T[j].f;
            i = R[f].b;
            uint32_t v = T[i].w;
            std::swap(T[j], T[i]);
            P[w] = i;
            P[v] = j;
            assert(T[i].f == f);
        }

        uint32_t& f = T[i].f;
        assert(f >= 2);
        if (R.find(f - 2) != R.cend()) {
            R[f - 2].e += 1;
        } else {
            R[f - 2] = {i, i + 1};  // one-past the end
        }
        R[f].b += 1;
        if (R[f].b == R[f].e) R.erase(f);
        f -= 2;
    }

    std::vector<wf> T;                         // array of pairs (weight,frequency)
    std::unordered_map<uint32_t, range> R;     // frequency --> [b,e) (range into T)
    std::unordered_map<uint32_t, uint32_t> P;  // weight --> p (position into T)
};

}  // namespace sshash