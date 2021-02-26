#include <string>
#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TCanvas.h"
#include <numeric>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include "include/coordinateTools.h"
#include "include/MultiCumulants/correlations/Types.hh"
#include "include/MultiCumulants/correlations/Result.hh"
#include "include/MultiCumulants/correlations/QVector.hh"
#include "include/MultiCumulants/correlations/recursive/FromQVector.hh"
#include "include/MultiCumulants/correlations/recurrence/FromQVector.hh"
#include "include/MultiCumulants/MultiCumulants/Subsets.h"
#include "include/MultiCumulants/MultiCumulants/QVector.h"
#include "include/MultiCumulants/MultiCumulants/QVectorSet.h"
#include "include/MultiCumulants/MultiCumulants/QTerms.h"
#include "include/MultiCumulants/MultiCumulants/Correlator.h"


void analyze( std::vector< std::string> files, std::string outputFileName = "output_cumulants_v2.root"){

  //stuff for qhat combination
  const int nqHats = 8;
  const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
  //units in pb
  float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
  TH1D * qHatHist = new TH1D("qhatHit",";;qHat",nqHats,qHatBoundaries);

  //vectors for loading branches
  float genQScale = 0;
  std::vector< float > * genJetPt = 0;
  std::vector< float > * genJetEta = 0;
  std::vector< float > * genJetPhi = 0;
  std::vector< int > * genJetChargedMultiplicity = 0;
  std::vector< std::vector< int > > * genDau_chg = 0;
  std::vector< std::vector< int > > * genDau_pid = 0;
  std::vector< std::vector< float > > * genDau_pt = 0;
  std::vector< std::vector< float > > * genDau_eta = 0;
  std::vector< std::vector< float > > * genDau_phi = 0;


  // Histograms
  TH1D*  hcN2 = new TH1D("hcN2",";N_{trk}",1000,0,1000);
  TH1D*  hcN4 = new TH1D("hcN4",";N_{trk}",1000,0,1000);
  TH1D*  hcN6 = new TH1D("hcN6",";N_{trk}",1000,0,1000);
  TH1D*  hcN8 = new TH1D("hcN8",";N_{trk}",1000,0,1000);

  //Init cumulants
  vector<double>  etasubmin_{-20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.};  //min eta of the tracks for subevents
  vector<double>  etasubmax_{20.,20.,20.,20.,20.,20.,20.,20.};  //max eta of the tracks for subevents
  vector<double>  ptsubmin_{0.,0.,0.,0.,0.,0.,0.,0.};  //min pt of the tracks for subevents
  vector<double>  ptsubmax_{3.,3.,3.,3.,3.,3.,3.,3.};  //max pt of the tracks for subevents  

  cumulant::Subset sub_1(2);
  sub_1.set(0, "pt", ptsubmin_[0], ptsubmax_[0]);
  sub_1.set(1, "eta", etasubmin_[0], etasubmax_[0]);
  cumulant::Subset sub_2(2);
  sub_2.set(0, "pt", ptsubmin_[1], ptsubmax_[1]);
  sub_2.set(1, "eta", etasubmin_[1], etasubmax_[1]);
  cumulant::Subset sub_3(2);
  sub_3.set(0, "pt", ptsubmin_[2], ptsubmax_[2]);
  sub_3.set(1, "eta", etasubmin_[2], etasubmax_[2]);
  cumulant::Subset sub_4(2);
  sub_4.set(0, "pt", ptsubmin_[3], ptsubmax_[3]);
  sub_4.set(1, "eta", etasubmin_[3], etasubmax_[3]);
  cumulant::Subset sub_5(2);
  sub_5.set(0, "pt", ptsubmin_[4], ptsubmax_[4]);
  sub_5.set(1, "eta", etasubmin_[4], etasubmax_[4]);
  cumulant::Subset sub_6(2);
  sub_6.set(0, "pt", ptsubmin_[5], ptsubmax_[5]);
  sub_6.set(1, "eta", etasubmin_[5], etasubmax_[5]);
  cumulant::Subset sub_7(2);
  sub_7.set(0, "pt", ptsubmin_[6], ptsubmax_[6]);
  sub_7.set(1, "eta", etasubmin_[6], etasubmax_[6]);
  cumulant::Subset sub_8(2);
  sub_8.set(0, "pt", ptsubmin_[7], ptsubmax_[7]);
  sub_8.set(1, "eta", etasubmin_[7], etasubmax_[7]);

  //Init sub-event method
  std::vector<int> harm_{2,2,2,2};     //harmonic order

  cumulant::Set set(8);
  set.setSubsetParams(0, sub_1);
  set.setSubsetParams(1, sub_2);
  set.setSubsetParams(2, sub_3);
  set.setSubsetParams(3, sub_4);
  set.setSubsetParams(4, sub_5);
  set.setSubsetParams(5, sub_6);
  set.setSubsetParams(6, sub_7);
  set.setSubsetParams(7, sub_8);

  correlations::HarmonicVector h(8);
  h[0] =  1*harm_[0];
  h[1] =  1*harm_[1];
  h[2] =  1*harm_[2];
  h[3] =  1*harm_[3];
  h[4] = -1*harm_[0];
  h[5] = -1*harm_[1];
  h[6] = -1*harm_[2];
  h[7] = -1*harm_[3];

  bool cweight_ = true;
  cumulant::QVectorSet qN_ = cumulant::QVectorSet(h, set, cweight_);

  double sumCN2_17_[1000]={0};
  double sumwCN2_17_[1000]={0};
  double sumCN4_51_[1000]={0};
  double sumwCN4_51_[1000]={0};
  double sumCN6_119_[1000]={0};
  double sumwCN6_119_[1000]={0};
  double sumCN8_[1000]={0};
  double sumwCN8_[1000]={0};

  // load phi efficiency correction map
  TH2D* heff_;
  if(cweight_)
  {
    TFile* feff_ = new TFile("output_phi.root");
    heff_ = (TH2D*)feff_->Get("hPhi");
    heff_->Scale(2*3.1416/heff_->Integral("width"));
  }

  std::cout << "calculating weighting factors for qhat bins" << std::endl;
  for(unsigned int f = 0; f<files.size(); f++){
    //for testing
    //if(f%10!=0) continue;
    //if(f%10==0) std::cout << f << "/" << files.size() << std::endl;
    std::cout << f << "/" << files.size() << std::endl;

    TFile * inputFile = TFile::Open(files.at(f).c_str(),"read");
    TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");

    t->SetBranchAddress("genQScale",&genQScale);
    t->GetEntry(0);

    qHatHist->Fill(genQScale,t->GetEntries());    

    inputFile->Close();
  } 
  qHatHist->Print("All");
  
  //analysis loop (loop over files first)
  for(unsigned int f = 0; f<files.size(); f++){
    //for testing
    //if(f%10!=0) continue;
    //if(f%10==0) std::cout << f << "/" << files.size() << std::endl;
    std::cout << f << "/" << files.size() << std::endl;

    //load file and set branches
    TFile * inputFile = TFile::Open(files.at(f).c_str(),"read");
    TTree * t = (TTree*)inputFile->Get("analyzer/trackTree");

    t->SetBranchAddress("genQScale",&genQScale);
    t->SetBranchAddress("genJetPt",&genJetPt);
    t->SetBranchAddress("genJetEta",&genJetEta);
    t->SetBranchAddress("genJetPhi",&genJetPhi);
    t->SetBranchAddress("genJetChargedMultiplicity",&genJetChargedMultiplicity); 
    t->SetBranchAddress("genDau_chg",&genDau_chg); 
    t->SetBranchAddress("genDau_pid",&genDau_pid); 
    t->SetBranchAddress("genDau_pt",&genDau_pt); 
    t->SetBranchAddress("genDau_eta",&genDau_eta); 
    t->SetBranchAddress("genDau_phi",&genDau_phi); 

    //event loop
    for(int i = 0; i<t->GetEntries(); i++){
      t->GetEntry(i);
      //weight by xsection/total number of gen events in the pthat bin
      float eventWeight = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));

      //jets loop
      for(int j = 0; j<genJetPt->size(); j++){
        //only take jets >500 GeV
        if(genJetPt->at(j) > 500){

          //constituent loop
          qN_.reset();
          std::vector<double> val(2,0.);
          int mult_ = 0;


          for(int k = 0; k<(genDau_chg->at(j)).size(); k++){
            //skip neutral particles and calculate variables with respect to jet axis:
            if( (genDau_chg->at(j)).at(k) == 0 ) continue;
            double ptStar = ptWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            double etaStar = etaWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            double phiStar = phiWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            double weight = 1;
            if(heff_) weight = 1./heff_->GetBinContent(heff_->FindBin(phiStar,genJetChargedMultiplicity->at(j)));

            //
            //insert constituent-level analysis code here
            //
            val[0] = ptStar;
            val[1] = etaStar;
            qN_.fill(val, phiStar, weight);
            if(ptStar<3 && ptStar>0 && fabs(etaStar)<20) mult_++;
          } //end of constituent loop

          // Calculate cumulants
          cumulant::QVectorMap& qNmap = qN_.getQ();

          cumulant::Correlator c2_17 = cumulant::Correlator(17, qNmap);
          double CN2_17_  = c2_17.v.real(); 
          double wCN2_17_ = c2_17.w;
          cumulant::Correlator c4_51 = cumulant::Correlator(51, qNmap);
          double CN4_51_  = c4_51.v.real();
          double wCN4_51_ = c4_51.w;
          cumulant::Correlator c6_119 = cumulant::Correlator(119, qNmap);
          double CN6_119_  = c6_119.v.real();
          double wCN6_119_ = c6_119.w;
          cumulant::Correlator c8 = cumulant::Correlator(255, qNmap);
          double CN8_  = c8.v.real();
          double wCN8_ = c8.w;

          if(wCN2_17_)  CN2_17_  = CN2_17_  / wCN2_17_;
          if(wCN4_51_)  CN4_51_  = CN4_51_  / wCN4_51_;
          if(wCN6_119_) CN6_119_ = CN6_119_ / wCN6_119_;
          if(wCN8_)     CN8_     = CN8_     / wCN8_;

          double cum2 = CN2_17_;
          double cum4 = CN4_51_ - 2*CN2_17_*CN2_17_;
          double cum6 = CN6_119_ - 9*CN4_51_*CN2_17_ + 12*CN2_17_*CN2_17_*CN2_17_;
          double cum8 = CN8_ - 16*CN6_119_*CN2_17_ - 18*CN4_51_*CN4_51_ + 144*CN4_51_*CN2_17_*CN2_17_ - 144*CN2_17_*CN2_17_*CN2_17_*CN2_17_;

          sumCN2_17_[mult_] += cum2*wCN2_17_;
          sumwCN2_17_[mult_] += wCN2_17_;
          sumCN4_51_[mult_] += cum4*wCN4_51_;
          sumwCN4_51_[mult_] += wCN4_51_;
          sumCN6_119_[mult_] += cum6*wCN6_119_;
          sumwCN6_119_[mult_] += wCN6_119_;
          sumCN8_[mult_] += cum8*wCN8_;
          sumwCN8_[mult_] += wCN8_;
        }
      }  //end of jet loop
    }  //end of event loop
    inputFile->Close();
  } //end of file loop 

  for(int i=0;i<hcN2->GetNbinsX();i++)
  {
    if(sumwCN2_17_[i])  sumCN2_17_[i]  = sumCN2_17_[i]  / sumwCN2_17_[i];
    if(sumwCN4_51_[i])  sumCN4_51_[i]  = sumCN4_51_[i]  / sumwCN4_51_[i];
    if(sumwCN6_119_[i]) sumCN6_119_[i] = sumCN6_119_[i] / sumwCN6_119_[i];
    if(sumwCN8_[i])     sumCN8_[i]     = sumCN8_[i]     / sumwCN8_[i];

    if(i>2) hcN2->SetBinContent(i+1,sumCN2_17_[i]);
    if(i>4) hcN4->SetBinContent(i+1,sumCN4_51_[i]);
    if(i>6) hcN6->SetBinContent(i+1,sumCN6_119_[i]);
    if(i>8) hcN8->SetBinContent(i+1,sumCN8_[i]);
  }

  TFile * outputFile = TFile::Open(outputFileName.data(),"recreate");
  hcN2->Write();
  hcN4->Write();
  hcN6->Write();
  hcN8->Write();
  outputFile->Close();
}

//Code enters execution here
int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList> <outputfile>" << std::endl;
    return 1;
  }  


  //read input parameters
  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());
  std::string outputFileName = argv[2];

  //read the file list and spit it into a vector of strings based on how the parallelization is to be done
  //each vector is a separate subset of the fileList based on the job number
  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      line++;
    }
  }

  analyze(listOfFiles,outputFileName);

  return 0; 
}
