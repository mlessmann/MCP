#include "barrier.h"

barrier::barrier(const std::size_t count) : m_count(count) {}

void barrier::wait() {
    std::unique_lock<std::mutex> lock(m_mutex);
    if (--m_count == 0)
        m_cv.notify_all();
    else
        m_cv.wait(lock, [this] {return m_count == 0;});
}
