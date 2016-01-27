#pragma once

#include <mutex>
#include <condition_variable>

class barrier {
public:
    barrier(const std::size_t count);
    void wait();

private:
    std::size_t             m_count;
    std::mutex              m_mutex;
    std::condition_variable m_cv;
};
