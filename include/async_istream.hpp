#pragma once

#include <condition_variable>
#include <deque>
#include <istream>
#include <mutex>
#include <thread>
#include <vector>

namespace sshash {

/*
    Read-ahead decorator for an std::istream: a producer thread reads (and,
    when the source is a decompressing stream, thereby decompresses) the
    source into a bounded queue of chunks, which this stream then serves.
    Decompression and parsing of the input then overlap with the queries
    instead of serializing with them. The source stream must outlive this
    object.
*/
struct async_istream : private std::streambuf, public std::istream  //
{
    static constexpr uint64_t chunk_bytes = 1 << 20;
    static constexpr uint64_t max_queued_chunks = 4;

    explicit async_istream(std::istream& source)
        : std::istream(this), m_stop(false), m_eof(false)  //
    {
        m_producer = std::thread([this, &source]() {
            while (true) {
                std::vector<char> chunk(chunk_bytes);
                source.read(chunk.data(), chunk.size());
                const uint64_t n = source.gcount();
                chunk.resize(n);
                std::unique_lock<std::mutex> lock(m_mutex);
                m_not_full.wait(lock,
                                [this] { return m_stop or m_queue.size() < max_queued_chunks; });
                if (m_stop) return;
                const bool done = n == 0 or !source;
                if (n != 0) m_queue.push_back(std::move(chunk));
                if (done) {
                    m_eof = true;
                    m_not_empty.notify_one();
                    return;
                }
                m_not_empty.notify_one();
            }
        });
    }

    ~async_istream() {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_stop = true;
        }
        m_not_full.notify_one();
        m_producer.join();
    }

private:
    int underflow() override {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_not_empty.wait(lock, [this] { return m_eof or !m_queue.empty(); });
        if (m_queue.empty()) return std::streambuf::traits_type::eof();
        m_current = std::move(m_queue.front());
        m_queue.pop_front();
        lock.unlock();
        m_not_full.notify_one();
        setg(m_current.data(), m_current.data(), m_current.data() + m_current.size());
        return std::streambuf::traits_type::to_int_type(*gptr());
    }

    std::thread m_producer;
    std::mutex m_mutex;
    std::condition_variable m_not_empty, m_not_full;
    std::deque<std::vector<char>> m_queue;
    std::vector<char> m_current;
    bool m_stop, m_eof;
};

}  // namespace sshash
