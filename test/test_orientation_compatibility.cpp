#include <filesystem>
#include <iterator>
#include <random>
#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"
using namespace sshash;
std::string rc(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) c = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'G' ? 'C' : 'G';
    return s;
}
template <class K, class O>
void check(std::string dir, std::string legacy, uint64_t k, std::string tag) {
    std::mt19937_64 rng(20261104);
    std::vector<std::string> words;
    for (int i = 0; i < 128; ++i) {
        std::string s;
        for (uint64_t j = 0; j < k; ++j) s += "ACGT"[rng() % 4];
        words.push_back(s);
    }
    auto input = dir + "/" + tag + ".fa";
    {
        std::ofstream out(input);
        for (auto s : words) out << ">x\n" << s << "\n";
    }
    build_configuration cfg;
    cfg.k = k;
    cfg.m = 15;
    cfg.verbose = false;
    cfg.tmp_dirname = dir;
    using D = dictionary<K, O>;
    D d;
    d.build(input, cfg);
    auto verify = [&](const D& x) {
        for (uint64_t i = 0; i < words.size(); ++i) {
            auto f = x.lookup(words[i].data(), false);
            auto r = x.lookup(rc(words[i]).data(), true);
            auto missing = x.lookup(rc(words[i]).data(), false);
            if (f.kmer_id != i || r.kmer_id != i || missing.kmer_id != constants::invalid_uint64)
                throw std::runtime_error("default or legacy mismatch");
        }
    };
    verify(d);
    auto check_version = [](const std::string& filename, uint8_t expected) {
        std::ifstream in(filename, std::ios::binary);
        essentials::generic_loader reader(in);
        essentials::version_number version(0, 0, 0);
        reader.visit(version);
        if (!in || version.x != expected) throw std::runtime_error("wrong serialized version");
    };
    auto read_bytes = [](const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), {});
    };
    auto path = dir + "/" + tag + ".sshash";
    essentials::save(d, path.c_str());
    check_version(path, constants::current_version_number::x);
    D loaded, mapped;
    essentials::load(loaded, path.c_str());
    essentials::mmap(mapped, path.c_str());
    verify(loaded);
    verify(mapped);
    if (!legacy.empty()) {
        auto oldpath = legacy + "/" + tag + ".sshash";
        check_version(oldpath, 6);
        essentials::load(loaded, oldpath.c_str());
        D oldmap;
        essentials::mmap(oldmap, oldpath.c_str());
        verify(loaded);
        verify(oldmap);
        for (auto const* source : {&loaded, &oldmap}) {
            auto upgraded = dir + "/" + tag + "-upgraded.sshash";
            essentials::save(*source, upgraded.c_str());
            check_version(upgraded, constants::current_version_number::x);
            if (read_bytes(upgraded) != read_bytes(oldpath))
                throw std::runtime_error("resaved v6 index differs from the original");
            D reloaded, remapped;
            essentials::load(reloaded, upgraded.c_str());
            essentials::mmap(remapped, upgraded.c_str());
            verify(reloaded);
            verify(remapped);
        }

        cfg.preserve_orientation = true;
        D oriented;
        oriented.build(input, cfg);
        auto np = dir + "/" + tag + "-oriented.sshash";
        essentials::save(oriented, np.c_str());
        check_version(np, constants::oriented_version_major);
        essentials::load(loaded, np.c_str());
        verify(loaded);
        auto resaved = dir + "/" + tag + "-oriented-resaved.sshash";
        essentials::save(loaded, resaved.c_str());
        if (read_bytes(resaved) != read_bytes(np))
            throw std::runtime_error("resaved oriented index differs from the original");
        essentials::load(loaded, oldpath.c_str());
        verify(loaded);
    }
}
int main(int argc, char** argv) {
    if (argc < 2) return 2;
    for (uint8_t major : {uint8_t(5), uint8_t(8)}) {
        bool rejected = false;
        try { util::check_version_number(essentials::version_number(major, 0, 0)); }
        catch (std::runtime_error const&) { rejected = true; }
        if (!rejected) throw std::runtime_error("unsupported version accepted");
    }
    std::filesystem::create_directory(argv[1]);
    std::string old = argc > 2 ? argv[2] : "";
#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
    if (!old.empty()) old += "/traditional";
#endif
    check<dna_uint_kmer_t<uint64_t>, decoded_offsets>(argv[1], old, 31, "k31-decoded");
    check<dna_uint_kmer_t<uint64_t>, encoded_offsets>(argv[1], old, 31, "k31-encoded");
    check<dna_uint_kmer_t<__uint128_t>, decoded_offsets>(argv[1], old, 41, "k41-decoded");
    check<dna_uint_kmer_t<__uint128_t>, encoded_offsets>(argv[1], old, 41, "k41-encoded");
    std::cout << "default/serialization checks passed\n";
}
