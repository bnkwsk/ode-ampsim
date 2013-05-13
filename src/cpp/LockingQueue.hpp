#pragma once

#include <condition_variable>
#include <mutex>
#include <queue>

template<typename T>
class LockingQueue
{
private:
    std::queue<T> queue;
    std::mutex mutex;
    std::condition_variable notEmptyCondition;
    std::condition_variable notFullCondition;
    const int maxSize = 100;
public:
    void push(const T &data)
    {
        std::unique_lock<std::mutex> lock(mutex);
        notFullCondition.wait(lock, [this] {return queue.size() != maxSize;});
        queue.push(data);
        lock.unlock();
        notEmptyCondition.notify_one();
    }
    T pop()
    {
        std::unique_lock<std::mutex> lock(mutex);
        T element;
        notEmptyCondition.wait(lock, [this] {return queue.size() != 0;});
        element = queue.front();
        queue.pop();
        lock.unlock();
        notFullCondition.notify_one();
        return element;
    }
};