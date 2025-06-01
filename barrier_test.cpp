#include <iostream>
#include <barrier>
#include <thread>

std::barrier sync_point(2);

void task(const std::string& name) {
    std::cout << name << " reached barrier\n";
    sync_point.arrive_and_wait();
    std::cout << name << " passed barrier\n";
}

int main() {
    std::thread t1(task, "Thread 1");
    std::thread t2(task, "Thread 2");
    t1.join();
    t2.join();
}
