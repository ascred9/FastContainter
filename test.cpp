#include <random>
#include <algorithm>
#include <chrono>

#include "FastContainer.h"

#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"

int main()
{
    // create randomer
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> udist(-200.0, 200.0);

    TGraph* gr_std = new TGraph(); 
    gr_std->SetName("gr_std");
    gr_std->SetTitle("Vector");
    gr_std->SetLineColor(kRed);
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

        // fill fast container
        auto startCreation = std::chrono::high_resolution_clock::now();
        FastContainer fc(vec, -200, 200);
        auto stopCreation = std::chrono::high_resolution_clock::now();
        auto durationCreation = std::chrono::duration_cast<std::chrono::microseconds>(stopCreation - startCreation);
        std::cout << "Creation time: " << durationCreation.count() << ", muSec" << std::endl;
    
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

        gr_std->AddPoint(N, durationStd.count());
        gr_fast->AddPoint(N, durationF.count());
        gr_fast_with_init->AddPoint(N, durationF.count()+durationCreation.count());
        gr_dens_std->AddPoint(N, durationStd.count() * 1./testN);
        gr_dens_fast->AddPoint(N, durationF.count() * 1./testN);
        gr_dens_fast_with_init->AddPoint(N, (durationF.count()+durationCreation.count()) * 1./testN);

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
    }

    TCanvas* c = new TCanvas("c", "Comparison", 1800, 900);
    c->Divide(2);

    c->cd(1)->SetGrid();
    c->cd(1)->SetLogy();
    c->cd(1);
    TMultiGraph* mg = new TMultiGraph("mg", "Comparison");
    mg->Add(gr_std);
    mg->Add(gr_fast);
    mg->Add(gr_fast_with_init);
    mg->GetXaxis()->SetTitle("Number of values");
    mg->GetYaxis()->SetTitle("Time, #muS");
    mg->Draw("AL");

    TLegend* legend = new TLegend(0.7, 0.6, 0.95, 0.7);
    legend->AddEntry("gr_std");
    legend->AddEntry("gr_fast");
    legend->AddEntry("gr_fast_with_init");
    legend->Draw();

    c->cd(2)->SetGrid();
    c->cd(2);
    TMultiGraph* mg_dens = new TMultiGraph("mg_dens", "Comparison");
    mg_dens->Add(gr_dens_std);
    mg_dens->Add(gr_dens_fast);
    mg_dens->Add(gr_dens_fast_with_init);
    mg_dens->GetXaxis()->SetTitle("Number of values");
    mg_dens->GetYaxis()->SetTitle("Time / Number of tests, #muS");
    mg_dens->Draw("AL");

    c->SaveAs("test.png");

    return 0;
}
