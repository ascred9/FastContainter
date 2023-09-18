#include <random>
#include <algorithm>
#include <chrono>

#include "FastContainer.h"

#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"

void testNearest(bool verbose)
{
    // create randomer
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> udist(-200.0, 200.0);

    TGraph* gr_std = new TGraph(); 
    gr_std->SetName("gr_std");
    gr_std->SetTitle("Vector");
    gr_std->SetLineColor(kRed);
    TGraph* gr_set = new TGraph(); 
    gr_set->SetName("gr_set");
    gr_set->SetTitle("Set");
    gr_set->SetLineColor(kGreen);
    TGraph* gr_fast = new TGraph(); 
    gr_fast->SetName("gr_fast");
    gr_fast->SetTitle("FastContainer");
    gr_fast->SetLineColor(kBlue);
    TGraph* gr_fast_with_init = new TGraph(); 
    gr_fast_with_init->SetName("gr_fast_with_init");
    gr_fast_with_init->SetTitle("FastContainerWithInit");
    gr_fast_with_init->SetLineColor(kBlack);

    TGraph* gr_dens_std = new TGraph(); 
    gr_dens_std->SetName("gr_dens_std");
    gr_dens_std->SetTitle("Vector");
    gr_dens_std->SetLineColor(kRed);
    TGraph* gr_dens_set = new TGraph(); 
    gr_dens_set->SetName("gr_dens_set");
    gr_dens_set->SetTitle("Set");
    gr_dens_set->SetLineColor(kGreen);
    TGraph* gr_dens_fast = new TGraph(); 
    gr_dens_fast->SetName("gr_dens_fast");
    gr_dens_fast->SetTitle("FastContainer");
    gr_dens_fast->SetLineColor(kBlue);
    TGraph* gr_dens_fast_with_init = new TGraph(); 
    gr_dens_fast_with_init->SetName("gr_dens_fast_with_init");
    gr_dens_fast_with_init->SetTitle("FastContainerWithInit");
    gr_dens_fast_with_init->SetLineColor(kBlack);

    int max_pow = 11;
    int testN = 1e5;

    for (int ipow = 0; ipow < max_pow; ++ipow)
    {
        // generate and fill input data
        int N = pow(2, ipow);
        std::cout << "Input number: " << N << std::endl;
        std::vector<std::pair<int, double>> vec;
        vec.reserve(N);
        for (int i=0; i<N; i++)
            vec.emplace_back(i, udist(gen));

        // generate test numbers
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
        if (verbose) std::cout << "Standart duration: " << durationStd.count() << ", muSec" << std::endl;

        auto startSet = std::chrono::high_resolution_clock::now();
        auto comp = [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs)
        {
            return lhs.second < rhs.second;
        };
        std::set<std::pair<int, double>, decltype(comp)> tmp_set (vec.begin(), vec.end(), comp);

        std::vector<int> resSet;
        resSet.reserve(testN);
        for (const auto& elem: test)
        {
            auto it = tmp_set.lower_bound({0, elem});
            auto pit = std::prev(it);
            double res = it->second - elem < elem - pit->second ? it->first : pit->first;
            resSet.push_back(res);
        }
        auto stopSet = std::chrono::high_resolution_clock::now();
        auto durationSet = std::chrono::duration_cast<std::chrono::microseconds>(stopSet - startSet);
        if (verbose) std::cout << "Standart duration: " << durationSet.count() << ", muSec" << std::endl;

        // TEST NEW SOLUTION
        // fill fast container
        auto startCreation = std::chrono::high_resolution_clock::now();
        FastContainer fc(-200, 200);
        fc.set(vec);

        std::vector<int> resF;
        resF.reserve(testN);
        auto startF = std::chrono::high_resolution_clock::now();
        for (const auto& elem: test)
        {
            int idx = *(fc.getClosestId(elem).first);
            resF.push_back(idx);
        }
        auto stopF = std::chrono::high_resolution_clock::now();
        auto durationF = std::chrono::duration_cast<std::chrono::microseconds>(stopF - startF);
        auto stopCreation = std::chrono::high_resolution_clock::now();
        auto durationCreation = std::chrono::duration_cast<std::chrono::microseconds>(stopCreation - startCreation);
        if (verbose){
            std::cout << "New duration: " << durationF.count() << ", muSec" << std::endl;
            std::cout << "With creation time: " << durationCreation.count() << ", muSec" << std::endl;
        }

        gr_std->AddPoint(N, durationStd.count());
        gr_set->AddPoint(N, durationSet.count());
        gr_fast->AddPoint(N, durationF.count());
        gr_fast_with_init->AddPoint(N, durationCreation.count());
        gr_dens_std->AddPoint(N, durationStd.count() * 1./testN);
        gr_dens_set->AddPoint(N, durationSet.count() * 1./testN);
        gr_dens_fast->AddPoint(N, durationF.count() * 1./testN);
        gr_dens_fast_with_init->AddPoint(N, durationCreation.count() * 1./testN);

        if (verbose){
            std::cout << "Ratio std / new: " << durationStd.count() * 1. / durationF.count() << std::endl;
            std::cout << "Ratio std / (new + cre): " << durationStd.count() * 1. / (durationF.count() + durationCreation.count()) << std::endl;
        }

        // compare values
        if (verbose) std::cout << "Check solutions" << std::endl;
        if (resStd.size() != resF.size())
            std::cout << "Different sizes" << std::endl;
        for (int i = 0; i < testN; i++)
        {
            if (resStd.at(i) == resF.at(i))
                continue;
    
            std::cout << test.at(i) << " \t" << resStd.at(i) << " " << vec.at(resStd.at(i)).second  << "\t" << test.at(i) - vec.at(resStd.at(i)).second << std::endl;
            std::cout << "\t\t" << resF.at(i) << " " << vec.at(resF.at(i)).second << "\t" << test.at(i) - vec.at(resF.at(i)).second << std::endl;
            std::cout << *(fc.getClosestId(test.at(i)).first) << std::endl;
        }
    }

    TCanvas* c = new TCanvas("c", "Comparison", 1800, 900);
    c->Divide(2);

    c->cd(1)->SetGrid();
    c->cd(1)->SetLogy();
    c->cd(1);
    TMultiGraph* mg = new TMultiGraph("mg", "Comparison getNearest");
    mg->Add(gr_std);
    mg->Add(gr_set);
    mg->Add(gr_fast);
    mg->Add(gr_fast_with_init);
    mg->GetXaxis()->SetTitle("Number of values");
    mg->GetYaxis()->SetTitle("Time, #muS");
    mg->Draw("AL");

    TLegend* legend = new TLegend(0.7, 0.6, 0.95, 0.7);
    legend->AddEntry("gr_std");
    legend->AddEntry("gr_set");
    legend->AddEntry("gr_fast");
    legend->AddEntry("gr_fast_with_init");
    legend->Draw();

    c->cd(2)->SetGrid();
    c->cd(2);
    TMultiGraph* mg_dens = new TMultiGraph("mg_dens", "Comparison getNearest, average time per step");
    mg_dens->Add(gr_dens_std);
    mg_dens->Add(gr_dens_set);
    mg_dens->Add(gr_dens_fast);
    mg_dens->Add(gr_dens_fast_with_init);
    mg_dens->GetXaxis()->SetTitle("Number of values");
    mg_dens->GetYaxis()->SetTitle("Time / Number of tests, #muS");
    mg_dens->Draw("AL");

    c->SaveAs("test.png");
}

void testRanges(bool verbose)
{
    // create randomer
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> udist(-200.0, 200.0);

    TGraph* gr_std = new TGraph(); 
    gr_std->SetName("gr_std");
    gr_std->SetTitle("Vector");
    gr_std->SetLineColor(kRed);
    TGraph* gr_set = new TGraph(); 
    gr_set->SetName("gr_multiset");
    gr_set->SetTitle("Multiset");
    gr_set->SetLineColor(kGreen);
    TGraph* gr_fast = new TGraph(); 
    gr_fast->SetName("gr_fast");
    gr_fast->SetTitle("FastContainer");
    gr_fast->SetLineColor(kBlue);
    TGraph* gr_fast_with_init = new TGraph(); 
    gr_fast_with_init->SetName("gr_fast_with_init");
    gr_fast_with_init->SetTitle("FastContainerWithInit");
    gr_fast_with_init->SetLineColor(kBlack);

    TGraph* gr_dens_std = new TGraph(); 
    gr_dens_std->SetName("gr_dens_std");
    gr_dens_std->SetTitle("Vector");
    gr_dens_std->SetLineColor(kRed);
    TGraph* gr_dens_set = new TGraph(); 
    gr_dens_set->SetName("gr_dens_multiset");
    gr_dens_set->SetTitle("Multiset");
    gr_dens_set->SetLineColor(kGreen);
    TGraph* gr_dens_fast = new TGraph(); 
    gr_dens_fast->SetName("gr_dens_fast");
    gr_dens_fast->SetTitle("FastContainer");
    gr_dens_fast->SetLineColor(kBlue);
    TGraph* gr_dens_fast_with_init = new TGraph(); 
    gr_dens_fast_with_init->SetName("gr_dens_fast_with_init");
    gr_dens_fast_with_init->SetTitle("FastContainerWithInit");
    gr_dens_fast_with_init->SetLineColor(kBlack);

    int max_pow = 11;
    int testN = 1e5;

    for (int ipow = 1; ipow < max_pow; ++ipow)
    {
        // generate and fill input data
        int N = pow(2, ipow);
        std::cout << "Input number: " << N << std::endl;
        std::vector<std::pair<int, double>> vec;
        vec.reserve(N);
        for (int i=0; i<N; i++)
            vec.emplace_back(i, udist(gen));

        // generate test numbers
        std::cout << "Test number: " << testN << std::endl;
        std::vector<double> test;
        test.reserve(testN);
        for (int i = 0; i<testN; i++)
        {
            test.push_back(udist(gen));
        }

        // TEST STD SOLUTION
        std::vector<std::vector<int>> resStd;
        resStd.reserve(testN);
        auto startStd = std::chrono::high_resolution_clock::now();
        for (int i=0; i<testN-1; ++i)
        {
            std::vector<int> tmp;
            double low = std::min(test[i], test[i+1]);
            double up = std::max(test[i], test[i+1]);
            std::for_each(vec.begin(), vec.end(), [&tmp, low, up](const auto& v){ if(v.second >= low && v.second <= up) tmp.push_back(v.first); });
            resStd.push_back(tmp);
        }
        auto stopStd = std::chrono::high_resolution_clock::now();
        auto durationStd = std::chrono::duration_cast<std::chrono::microseconds>(stopStd - startStd);
        if (verbose) std::cout << "Standart duration: " << durationStd.count() << ", muSec" << std::endl;

        auto startSet = std::chrono::high_resolution_clock::now();
        auto comp = [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs)
        {
            return lhs.second < rhs.second;
        };
        std::set<std::pair<int, double>, decltype(comp)> tmp_multiset (vec.begin(), vec.end(), comp);

        std::vector<std::vector<int>> resSet;
        resSet.reserve(testN);
        for (int i=0; i<testN-1; ++i)
        {
            // simple
            double low = std::min(test[i], test[i+1]);
            double up = std::max(test[i], test[i+1]);

            const auto lit = tmp_multiset.lower_bound({0, low});
            const auto uit = tmp_multiset.upper_bound({0, up});
            std::vector<int> tmp;
            tmp.reserve(std::distance(lit, uit)+1);
            std::for_each(lit, uit, [&tmp](const auto& v){ tmp.push_back(v.first); });
            
            resSet.push_back(tmp);
        }
        auto stopSet = std::chrono::high_resolution_clock::now();
        auto durationSet = std::chrono::duration_cast<std::chrono::microseconds>(stopSet - startSet);
        if (verbose) std::cout << "Standart duration: " << durationSet.count() << ", muSec" << std::endl;

        // TEST NEW SOLUTION
        // fill fast container
        auto startCreation = std::chrono::high_resolution_clock::now();
        FastContainer fc(-200, 200);
        fc.set(vec);

        std::vector<std::vector<int>> resF;
        resF.reserve(testN);
        auto startF = std::chrono::high_resolution_clock::now();
        for (int i=0; i<testN-1; ++i)
        {
            double low = std::min(test[i], test[i+1]);
            double up = std::max(test[i], test[i+1]);
            const auto& [fit, lit] = fc.getIdsInRange(low, up);
            auto tmp = std::vector<int>(fit, lit);
            resF.push_back(tmp);
        }
        auto stopF = std::chrono::high_resolution_clock::now();
        auto durationF = std::chrono::duration_cast<std::chrono::microseconds>(stopF - startF);
        auto stopCreation = std::chrono::high_resolution_clock::now();
        auto durationCreation = std::chrono::duration_cast<std::chrono::microseconds>(stopCreation - startCreation);
        if (verbose){
            std::cout << "New duration: " << durationF.count() << ", muSec" << std::endl;
            std::cout << "With creation time: " << durationCreation.count() << ", muSec" << std::endl;
        }

        gr_std->AddPoint(N, durationStd.count());
        gr_set->AddPoint(N, durationSet.count());
        gr_fast->AddPoint(N, durationF.count());
        gr_fast_with_init->AddPoint(N, durationCreation.count());
        gr_dens_std->AddPoint(N, durationStd.count() * 1./testN);
        gr_dens_set->AddPoint(N, durationSet.count() * 1./testN);
        gr_dens_fast->AddPoint(N, durationF.count() * 1./testN);
        gr_dens_fast_with_init->AddPoint(N, durationCreation.count() * 1./testN);

        if (verbose){
            std::cout << "Ratio std / new: " << durationStd.count() * 1. / durationF.count() << std::endl;
            std::cout << "Ratio std / (new + cre): " << durationStd.count() * 1. / (durationF.count() + durationCreation.count()) << std::endl;
        }

        // compare values
        if (verbose) std::cout << "Check solutions" << std::endl;
        if (resStd.size() != resF.size())
            std::cout << "Different sizes" << std::endl;
        for (int i = 0; i < testN-1; i++)
        {
            bool flag1 = true;
            bool flag2 = true;
            if (resStd.at(i).size() == resF.at(i).size())
            {
                flag1 = std::any_of(resStd.at(i).begin(), resStd.at(i).end(),
                    [&resF, i](const int v){return std::find(resF.at(i).begin(), resF.at(i).end(), v) == resF.at(i).end();});
                
                flag2 = std::any_of(resF.at(i).begin(), resF.at(i).end(),
                    [&resStd, i](const int v){return std::find(resStd.at(i).begin(), resStd.at(i).end(), v) == resStd.at(i).end();});

                if (!flag1 && !flag2)
                    continue;
                
                std::cout << i << "Error: " << flag1 << " " << flag2 << std::endl; 
            }
            else
                std::cout << i << " Different size: " << resStd.at(i).size() << " " << resF.at(i).size() << std::endl; 

            double low = std::min(test[i], test[i+1]);
            double up = std::max(test[i], test[i+1]);
            std::cout << "lowerZ: " << low << " upperZ: " << up << std::endl;
            std::cout << "Print input: ";
            std::for_each(vec.begin(), vec.end(), [](const auto& v){std::cout << v.first << "-" << v.second << " ";});
            std::cout << std::endl;
            std::cout << "Print std: ";
            std::for_each(resStd.at(i).begin(), resStd.at(i).end(), [](const auto& v){std::cout << v << " ";});
            std::cout << std::endl;
            
            std::cout << "Print fast: ";
            std::for_each(resF.at(i).begin(), resF.at(i).end(), [](const auto& v){std::cout << v << " ";});
            std::cout << std::endl;

            break;
        }
    }

    TCanvas* c = new TCanvas("c", "Comparison", 1800, 900);
    c->Divide(2);

    c->cd(1)->SetGrid();
    c->cd(1)->SetLogy();
    c->cd(1);
    TMultiGraph* mg = new TMultiGraph("mg", "Comparison getRange");
    mg->Add(gr_std);
    mg->Add(gr_set);
    mg->Add(gr_fast);
    mg->Add(gr_fast_with_init);
    mg->GetXaxis()->SetTitle("Number of values");
    mg->GetYaxis()->SetTitle("Time, #muS");
    mg->Draw("AL");

    TLegend* legend = new TLegend(0.7, 0.6, 0.95, 0.7);
    legend->AddEntry("gr_std");
    legend->AddEntry("gr_multiset");
    legend->AddEntry("gr_fast");
    legend->AddEntry("gr_fast_with_init");
    legend->Draw();

    c->cd(2)->SetGrid();
    c->cd(2);
    TMultiGraph* mg_dens = new TMultiGraph("mg_dens", "Comparison getRange, average time per step");
    mg_dens->Add(gr_dens_std);
    mg_dens->Add(gr_dens_set);
    mg_dens->Add(gr_dens_fast);
    mg_dens->Add(gr_dens_fast_with_init);
    mg_dens->GetXaxis()->SetTitle("Number of values");
    mg_dens->GetYaxis()->SetTitle("Time / Number of tests, #muS");
    mg_dens->Draw("AL");

    c->SaveAs("testRanges.png");
}

int main()
{
    testNearest(false);
    testRanges(false);
    return 0;
}
