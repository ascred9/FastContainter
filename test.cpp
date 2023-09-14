#include <random>
#include <algorithm>
#include <chrono>

#include "FastContainer.h"

int main()
{
    // create randomer
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> udist(-200.0, 200.0);

    // generate and fill input data
    int N = 1e2;
    std::cout << "Input number: " << N << std::endl;
    std::vector<std::pair<int, double>> vec;
    vec.reserve(N);
    for (int i=0; i<N; i++)
        vec.emplace_back(i, udist(gen));

    // fill fast container
    auto startCreation = std::chrono::high_resolution_clock::now();
    FastContainer fc(vec, -200, 200);
    auto stopCreation = std::chrono::high_resolution_clock::now();
    auto durationCreation = std::chrono::duration_cast<std::chrono::microseconds>(stopCreation - startCreation);
    std::cout << "Creation time: " << durationCreation.count() << ", muSec" << std::endl;
    
    // generate test numbers
    int testN = 1e5;
    std::cout << "Test number: " << testN << std::endl;
    std::vector<double> test;
    test.reserve(testN);
    for (int i = 0; i<testN; i++)
    {
        test.push_back(udist(gen));
    }

    // TEST STD SOLUTION
    std::vector<int> resStd;
    resStd.reserve(testN);
    auto startStd = std::chrono::high_resolution_clock::now();
    for (const auto& elem: test)
    {
        auto it = std::min_element(vec.begin(), vec.end(), [&elem](const auto& lhs, const auto& rhs){
            return std::abs(lhs.second-elem) < std::abs(rhs.second-elem);
        });
        resStd.push_back(it->first);
    }
    auto stopStd = std::chrono::high_resolution_clock::now();
    auto durationStd = std::chrono::duration_cast<std::chrono::microseconds>(stopStd - startStd);
    std::cout << "Standart duration: " << durationStd.count() << ", muSec" << std::endl;

    // TEST NEW SOLUTION
    std::vector<int> resF;
    resF.reserve(testN);
    auto startF = std::chrono::high_resolution_clock::now();
    for (const auto& elem: test)
    {
        int idx = fc.getClosestId(elem);
        resF.push_back(idx);
    }
    auto stopF = std::chrono::high_resolution_clock::now();
    auto durationF = std::chrono::duration_cast<std::chrono::microseconds>(stopF - startF);
    std::cout << "New duration: " << durationF.count() << ", muSec" << std::endl;

    std::cout << "Ratio std / new: " << durationStd.count() * 1. / durationF.count() << std::endl;
    std::cout << "Ratio std / (new + cre): " << durationStd.count() * 1. / (durationF.count() + durationCreation.count()) << std::endl;

    // compare values
    std::cout << "Check solutions" << std::endl;
    if (resStd.size() != resF.size())
        std::cout << "Different sizes" << std::endl;
    for (int i = 0; i < testN; i++)
    {
        if (resStd.at(i) == resF.at(i))
            continue;
    
        std::cout << test.at(i) << " \t" << resStd.at(i) << " " << vec.at(resStd.at(i)).second  << "\t" << test.at(i) - vec.at(resStd.at(i)).second << std::endl;
        std::cout << "\t\t" << resF.at(i) << " " << vec.at(resF.at(i)).second << "\t" << test.at(i) - vec.at(resF.at(i)).second << std::endl;
        std::cout << fc.getClosestId(test.at(i)) << std::endl;
    }

    return 0;
}
