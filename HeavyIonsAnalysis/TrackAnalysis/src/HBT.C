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

void analyze( std::vector< std::string> files){

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

  // Histograms for HBT analysis
  TH1D* hQSignal = new TH1D("hQSignal",";q_{inv}",1000,0,10);
  TH1D* hQBkg = new TH1D("hQBkg",";q_{inv}",1000,0,10);

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

    std::vector< std::vector< float > > ptStar_pos_vectvect;
    std::vector< std::vector< float > > etaStar_pos_vectvect;
    std::vector< std::vector< float > > phiStar_pos_vectvect;
    std::vector< std::vector< float > > ptStar_neg_vectvect;
    std::vector< std::vector< float > > etaStar_neg_vectvect;
    std::vector< std::vector< float > > phiStar_neg_vectvect;
    
    for(int i = 0; i<t->GetEntries(); i++){
      t->GetEntry(i);
      //weight by xsection/total number of gen events in the pthat bin
      float eventWeight = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));

      //jets loop
      for(int j = 0; j<genJetPt->size(); j++){
        //only take jets >500 GeV
        if(genJetPt->at(j) > 500){

          //constituent loop

          std::vector< float > ptStar_pos_vect;
          std::vector< float > etaStar_pos_vect;
          std::vector< float > phiStar_pos_vect;
          std::vector< float > ptStar_neg_vect;
          std::vector< float > etaStar_neg_vect;
          std::vector< float > phiStar_neg_vect;
          for(int k = 0; k<(genDau_chg->at(j)).size(); k++){
            //skip neutral particles and calculate variables with respect to jet axis:
            if( (genDau_chg->at(j)).at(k) == 0 ) continue;
            float ptStar = ptWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float etaStar = etaWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float phiStar = phiWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));

            if(ptStar<0 || ptStar>3) continue;
 
            if( (genDau_chg->at(j)).at(k)>0 )
            {
              ptStar_pos_vect.push_back(ptStar);
              etaStar_pos_vect.push_back(etaStar);
              phiStar_pos_vect.push_back(phiStar);
            }
            if( (genDau_chg->at(j)).at(k)<0 )
            {
              ptStar_neg_vect.push_back(ptStar);
              etaStar_neg_vect.push_back(etaStar);
              phiStar_neg_vect.push_back(phiStar);
            }
            //
            //insert constituent-level analysis code here
            //

          } //end of constituent loop
         
          // fill signal pairs
          for(int ii = 0; ii<(int)ptStar_pos_vect.size(); ii++)
            for(int jj = ii+1; jj<(int)ptStar_pos_vect.size(); jj++)
            {
              TVector3 vect1, vect2;
              vect1.SetPtEtaPhi(ptStar_pos_vect.at(ii),etaStar_pos_vect.at(ii),phiStar_pos_vect.at(ii));
              vect2.SetPtEtaPhi(ptStar_pos_vect.at(jj),etaStar_pos_vect.at(jj),phiStar_pos_vect.at(jj));
 
              hQSignal->Fill((vect1-vect2).Mag());
            }

          for(int ii = 0; ii<(int)ptStar_neg_vect.size(); ii++)
            for(int jj = ii+1; jj<(int)ptStar_neg_vect.size(); jj++)
            {
              TVector3 vect1, vect2;
              vect1.SetPtEtaPhi(ptStar_neg_vect.at(ii),etaStar_neg_vect.at(ii),phiStar_neg_vect.at(ii));
              vect2.SetPtEtaPhi(ptStar_neg_vect.at(jj),etaStar_neg_vect.at(jj),phiStar_neg_vect.at(jj));

              hQSignal->Fill((vect1-vect2).Mag());
            }

          // fill bkg pairs
          for(int mm = 0; mm<(int)ptStar_pos_vectvect.size(); mm++)
            for(int ii = 0; ii<(int)ptStar_pos_vect.size(); ii++)
              for(int jj = 0; jj<(int)(ptStar_pos_vectvect.at(mm)).size(); jj++)
              {
                TVector3 vect1, vect2;
                vect1.SetPtEtaPhi(ptStar_pos_vect.at(ii),etaStar_pos_vect.at(ii),phiStar_pos_vect.at(ii));
                vect2.SetPtEtaPhi((ptStar_pos_vectvect.at(mm)).at(jj),(etaStar_pos_vectvect.at(mm)).at(jj),(phiStar_pos_vectvect.at(mm)).at(jj));

                hQBkg->Fill((vect1-vect2).Mag());
              }          

          for(int mm = 0; mm<(int)ptStar_neg_vectvect.size(); mm++)
            for(int ii = 0; ii<(int)ptStar_neg_vect.size(); ii++)
              for(int jj = 0; jj<(int)(ptStar_neg_vectvect.at(mm)).size(); jj++)
              {
                TVector3 vect1, vect2;
                vect1.SetPtEtaPhi(ptStar_neg_vect.at(ii),etaStar_neg_vect.at(ii),phiStar_neg_vect.at(ii));
                vect2.SetPtEtaPhi((ptStar_neg_vectvect.at(mm)).at(jj),(etaStar_neg_vectvect.at(mm)).at(jj),(phiStar_neg_vectvect.at(mm)).at(jj));

                hQBkg->Fill((vect1-vect2).Mag());
              }

          // fill mixed event vectors
          if((int)ptStar_pos_vectvect.size()==10)
          {
            ptStar_pos_vectvect.erase(ptStar_pos_vectvect.begin());
            etaStar_pos_vectvect.erase(etaStar_pos_vectvect.begin());
            phiStar_pos_vectvect.erase(phiStar_pos_vectvect.begin());
          }
          ptStar_pos_vectvect.push_back(ptStar_pos_vect);
          etaStar_pos_vectvect.push_back(etaStar_pos_vect);
          phiStar_pos_vectvect.push_back(phiStar_pos_vect);

          if((int)ptStar_neg_vectvect.size()==10)
          {
            ptStar_neg_vectvect.erase(ptStar_neg_vectvect.begin());
            etaStar_neg_vectvect.erase(etaStar_neg_vectvect.begin());
            phiStar_neg_vectvect.erase(phiStar_neg_vectvect.begin());
          }
          ptStar_neg_vectvect.push_back(ptStar_neg_vect);
          etaStar_neg_vectvect.push_back(etaStar_neg_vect);
          phiStar_neg_vectvect.push_back(phiStar_neg_vect);

        }
      }  //end of jet loop
    }  //end of event loop
    inputFile->Close();
  } //end of file loop
 
  TFile * outputFile = TFile::Open("output_HBT.root","recreate");
  hQSignal->Write();
  hQBkg->Write();
  outputFile->Close();
}

//Code enters execution here
int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList>" << std::endl;
    return 1;
  }  


  //read input parameters
  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());


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

  analyze(listOfFiles);

  return 0; 
}
