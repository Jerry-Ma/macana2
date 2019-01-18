#include <benchmark/benchmark.h>
#include <gtest/gtest.h>
#include <logging.h>

int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    benchmark::Initialize(&argc, argv);
    std::cout << "Running tests:\n";

    int result = RUN_ALL_TESTS();
    if (!testing::FLAGS_gtest_list_tests) {
        std::cout << "\nRunning benchmarks:\n";
        spdlog::set_level(spdlog::level::critical);
        benchmark::RunSpecifiedBenchmarks();
    }
    return result;
}