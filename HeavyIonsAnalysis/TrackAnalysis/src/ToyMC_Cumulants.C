#include <string>
#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
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

#define PI 3.1416

void analyze(int nevents, std::string outputFileName = "output_toyMC.root"){

  // Histograms  
  TH1D*  hcN1 = new TH1D("hcN1",";N_{trk}",1000,0,1000);
  TH1D*  hcN2 = new TH1D("hcN2",";N_{trk}",1000,0,1000);
  TH1D*  hcN4 = new TH1D("hcN4",";N_{trk}",1000,0,1000);
  TH1D*  hcN6 = new TH1D("hcN6",";N_{trk}",1000,0,1000);
  TH1D*  hcN8 = new TH1D("hcN8",";N_{trk}",1000,0,1000);

  //Init cumulants
//  4-subevents, R-0.8
//  vector<double>  etasubmin_{2.,0.8,2.,0.8, 3. ,4.5,3. ,4.5};  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{3.,2. ,3.,2. , 4.5,7. ,4.5,7. };  //max eta of the tracks for subevents

//  3-subevents, R-0.8
//  vector<double>  etasubmin_{2.5,2.5,0.8,0.8,0.8,4. ,4. ,4. };  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{4., 4. ,2.5,2.5,2.5,7. ,7. ,7. };  //max eta of the tracks for subevents

//  2-subevents, R-0.8
//  vector<double>  etasubmin_{0.8,0.8,0.8,0.8,3. ,3. ,3. ,3. };  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{3. ,3. ,3. ,3. ,7. ,7. ,7. ,7. };  //max eta of the tracks for subevents

//  standard, R-0.8
//  vector<double>  etasubmin_{0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8};  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{7. ,7. ,7. ,7. ,7. ,7. ,7. ,7. };  //max eta of the tracks for subevents

//  4-subevents, R-0.4
//  vector<double>  etasubmin_{2. ,1.6,2.,1.6, 3. ,4.5,3. ,4.5};  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{3. ,2. ,3.,2. , 4.5,7. ,4.5,7. };  //max eta of the tracks for subevents

//  3-subevents, R-0.4
  vector<double>  etasubmin_{2.5,2.5,1.6,0.0,1.6,4. ,4. ,4. };  //min eta of the tracks for subevents
  vector<double>  etasubmax_{4., 4. ,2.5,0.0,2.5,7. ,7. ,7. };  //max eta of the tracks for subevents

//  2-subevents, R-0.4
//  vector<double>  etasubmin_{1.6,1.6,1.6,1.6,3. ,3. ,3. ,3. };  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{3. ,3. ,3. ,3. ,7. ,7. ,7. ,7. };  //max eta of the tracks for subevents

//  standard, R-0.4
//  vector<double>  etasubmin_{1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6};  //min eta of the tracks for subevents
//  vector<double>  etasubmax_{7. ,7. ,7. ,7. ,7. ,7. ,7. ,7. };  //max eta of the tracks for subevents

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

  bool cweight_ = false;
  cumulant::QVectorSet qN_ = cumulant::QVectorSet(h, set, cweight_);

  double sumCN1_1_[1000]={0};
  double sumwCN1_1_[1000]={0};
  double sumCN2_17_[1000]={0};
  double sumwCN2_17_[1000]={0};
  double sumCN2_18_[1000]={0};
  double sumwCN2_18_[1000]={0};
  double sumCN2_33_[1000]={0};
  double sumwCN2_33_[1000]={0};
  double sumCN2_34_[1000]={0};
  double sumwCN2_34_[1000]={0};
  double sumCN4_51_[1000]={0};
  double sumwCN4_51_[1000]={0};
  double sumCN6_119_[1000]={0};
  double sumwCN6_119_[1000]={0};
  double sumCN8_[1000]={0};
  double sumwCN8_[1000]={0};

  gRandom->SetSeed(0);
  TF1* func = new TF1("func","1+2*[0]*cos(2*(x-[1]))",-PI,PI);

  // load phi efficiency correction map
  TH2D* heff_;
  if(cweight_)
  {
    TFile* feff_ = new TFile("EfficiencyMap.root");
    heff_ = (TH2D*)feff_->Get("hEtaPhi");
//    heff_->Scale(2*3.1416/heff_->Integral("width"));
  }

  TFile* fin_ = new TFile("EfficiencyMap.root");
  TH2D* hJetPtMultTmp = (TH2D*)fin_->Get("hJetPtMult");
  TH1D* hMultTmp = (TH1D*)hJetPtMultTmp->ProjectionY("hMultTmp",-1,-1,"e"); 
  TH2D* hEtaPhiTmp = (TH2D*)fin_->Get("hEtaPhi");  

  //analysis loop
  for(unsigned int f = 0; f<nevents; f++){
    if(!(f % 10000)) std::cout << f << "/" << nevents << std::endl;

    qN_.reset();
    std::vector<double> val(2,0.);
    unsigned int mult_ = hMultTmp->GetRandom();
//    unsigned int mult_ = 20;

    func->SetParameter(0,0.3);
    func->SetParameter(1,2*(gRandom->Rndm()-0.5)*PI);

    for(unsigned int m = 0; m<mult_; m++){

            double ptStar = 0.5; 
            double etaStar;
            double phiStar; 
            hEtaPhiTmp->GetRandom2(phiStar,etaStar);
            phiStar = func->GetRandom();
            
            double weight = 1.;
            if(heff_) weight = weight / heff_->GetBinContent(heff_->FindBin(phiStar,etaStar));

            //
            //insert constituent-level analysis code here
            //
            val[0] = ptStar;
            val[1] = etaStar;
            qN_.fill(val, phiStar, weight);
    } //end of constituent loop
 
          // Calculate cumulants

          cumulant::QVectorMap& qNmap = qN_.getQ();

          cumulant::Correlator c1_1 = cumulant::Correlator(1, qNmap);
          double CN1_1_  = c1_1.v.real();
          double wCN1_1_ = c1_1.w;
          cumulant::Correlator c2_17 = cumulant::Correlator(17, qNmap);
          double CN2_17_  = c2_17.v.real(); 
          double wCN2_17_ = c2_17.w;
          cumulant::Correlator c2_18 = cumulant::Correlator(18, qNmap);
          double CN2_18_  = c2_18.v.real();
          double wCN2_18_ = c2_18.w;
          cumulant::Correlator c2_33 = cumulant::Correlator(33, qNmap);
          double CN2_33_  = c2_33.v.real();
          double wCN2_33_ = c2_33.w;
          cumulant::Correlator c2_34 = cumulant::Correlator(34, qNmap);
          double CN2_34_  = c2_34.v.real();
          double wCN2_34_ = c2_34.w;
          cumulant::Correlator c4_51 = cumulant::Correlator(51, qNmap);
          double CN4_51_  = c4_51.v.real();
          double wCN4_51_ = c4_51.w;
          cumulant::Correlator c6_119 = cumulant::Correlator(119, qNmap);
          double CN6_119_  = c6_119.v.real();
          double wCN6_119_ = c6_119.w;
          cumulant::Correlator c8 = cumulant::Correlator(255, qNmap);
          double CN8_  = c8.v.real();
          double wCN8_ = c8.w;

          sumCN1_1_[mult_] += CN1_1_;
          sumwCN1_1_[mult_] += wCN1_1_;
          sumCN2_17_[mult_] += CN2_17_;
          sumwCN2_17_[mult_] += wCN2_17_;
          sumCN2_18_[mult_] += CN2_18_;
          sumwCN2_18_[mult_] += wCN2_18_;
          sumCN2_33_[mult_] += CN2_33_;
          sumwCN2_33_[mult_] += wCN2_33_;
          sumCN2_34_[mult_] += CN2_34_;
          sumwCN2_34_[mult_] += wCN2_34_;
          sumCN4_51_[mult_] += CN4_51_;
          sumwCN4_51_[mult_] += wCN4_51_;
          sumCN6_119_[mult_] += CN6_119_;
          sumwCN6_119_[mult_] += wCN6_119_;
          sumCN8_[mult_] += CN8_;
          sumwCN8_[mult_] += wCN8_;
  } //end of event loop 

  for(int i=0;i<hcN2->GetNbinsX();i++)
  {
    if(sumwCN1_1_[i])  sumCN1_1_[i]  = sumCN1_1_[i]  / sumwCN1_1_[i];
    if(sumwCN2_17_[i])  sumCN2_17_[i]  = sumCN2_17_[i]  / sumwCN2_17_[i];
    if(sumwCN2_18_[i])  sumCN2_18_[i]  = sumCN2_18_[i]  / sumwCN2_18_[i];
    if(sumwCN2_33_[i])  sumCN2_33_[i]  = sumCN2_33_[i]  / sumwCN2_33_[i];
    if(sumwCN2_34_[i])  sumCN2_34_[i]  = sumCN2_34_[i]  / sumwCN2_34_[i];
    if(sumwCN4_51_[i])  sumCN4_51_[i]  = sumCN4_51_[i]  / sumwCN4_51_[i];
    if(sumwCN6_119_[i]) sumCN6_119_[i] = sumCN6_119_[i] / sumwCN6_119_[i];
    if(sumwCN8_[i])     sumCN8_[i]     = sumCN8_[i]     / sumwCN8_[i];

    double cum1 = sumCN1_1_[i];
    double cum2 = sumCN2_17_[i];
    double cum4 = sumCN4_51_[i] - sumCN2_17_[i]*sumCN2_34_[i] - sumCN2_18_[i]*sumCN2_33_[i];
    double cum6 = sumCN6_119_[i] - 9*sumCN4_51_[i]*sumCN2_17_[i] + 12*sumCN2_17_[i]*sumCN2_17_[i]*sumCN2_17_[i];
    double cum8 = sumCN8_[i] - 16*sumCN6_119_[i]*sumCN2_17_[i] - 18*sumCN4_51_[i]*sumCN4_51_[i] + 144*sumCN4_51_[i]*sumCN2_17_[i]*sumCN2_17_[i] - 144*sumCN2_17_[i]*sumCN2_17_[i]*sumCN2_17_[i]*sumCN2_17_[i];

    hcN1->SetBinContent(i+1,cum1);
    if(i>2) hcN2->SetBinContent(i+1,cum2);
    if(i>4) hcN4->SetBinContent(i+1,cum4);
    if(i>6) hcN6->SetBinContent(i+1,cum6);
    if(i>8) hcN8->SetBinContent(i+1,cum8);
  }

  TFile * outputFile = TFile::Open(outputFileName.data(),"recreate");
  hcN1->Write();
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
    std::cout << "Usage: Z_mumu_Channel <nevents> <outputfile>" << std::endl;
    return 1;
  }

  //read input parameters
  int nevents;
  sscanf (argv[1],"%d",&nevents);

  std::string outputFileName = argv[2];

  analyze(nevents,outputFileName);

  return 0; 
}
