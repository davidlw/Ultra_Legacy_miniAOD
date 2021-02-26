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

void analyze( std::vector< std::string> files, std::string outputFileName = "output_photons.root"){

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
  std::vector< std::vector< int > > * genDau_mom = 0;
  std::vector< std::vector< float > > * genDau_pt = 0;
  std::vector< std::vector< float > > * genDau_eta = 0;
  std::vector< std::vector< float > > * genDau_phi = 0;

  // Histograms for photon analysis
  TH1D* hPt = new TH1D("hPt",";p_{T}",1000,0,10);
  TH1D* hPtFrag = new TH1D("hPtFrag",";p_{T}",1000,0,10);
  TH2D* hPtBvsPtJFrag = new TH2D("hPtBvsPtJFrag",";p^{Jet}_{T};p^{Beam}_{T}",1000,0,10,2000,0,1000);
  TH2D* hEtaBvsEtaJFrag = new TH2D("hEtaBvsEtaJFrag",";#eta^{Jet}_{T};#eta^{Beam}_{T}",200,-10,10,200,-10,10);
  TH2D* hPtvsMult = new TH2D("hPtvsMult",";nMult;p_{T}",300,0,300,1000,0,10);
  TH2D* hPtFragvsMult = new TH2D("hPtFragvsMult",";nMult;p_{T}",300,0,300,1000,0,10);
  TH1D* hPID = new TH1D("hPID",";pdgcode",60000,-30000,30000);
  TH2D* hPIDvsMult = new TH2D("hPIDvsMult",";nMult;pdgcode",300,0,300,60000,-30000,30000);
  TH1D* hMomPID = new TH1D("hMomPID",";mother pdgcode",10000,-5000,5000);
  TH2D* hMomPIDvsMult = new TH2D("hMomPIDvsMult",";nMult;mother pdgcode",300,0,300,10000,-5000,5000);
  TH2D* hJetPvsMult = new TH2D("hJetPvsMult",";nMult;p (GeV)",300,0,300,200,0,2000);


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
    t->SetBranchAddress("genDau_mom",&genDau_mom);

    //event loop
    for(int i = 0; i<t->GetEntries(); i++){
      t->GetEntry(i);
      //weight by xsection/total number of gen events in the pthat bin
      float eventWeight = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));

      //jets loop
      for(int j = 0; j<genJetPt->size(); j++){
        //only take jets >500 GeV
        if(genJetPt->at(j) > 500 && fabs(genJetEta->at(j)) < 1.6){

          int nmult = genJetChargedMultiplicity->at(j);
          hJetPvsMult->Fill(nmult, genJetPt->at(j)*cosh(genJetEta->at(j)));

          //constituent loop
          for(int k = 0; k<nmult; k++){
            //skip neutral particles and calculate variables with respect to jet axis:
//            if( (genDau_chg->at(j)).at(k) != 0 ) continue;
            float pt = (genDau_pt->at(j)).at(k);
            float eta = (genDau_eta->at(j)).at(k);
            float ptStar = ptWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float etaStar = etaWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));
            float phiStar = phiWRTJet(genJetPt->at(j), genJetEta->at(j), genJetPhi->at(j), (genDau_pt->at(j)).at(k), (genDau_eta->at(j)).at(k), (genDau_phi->at(j)).at(k));

        //    if(ptStar<0 || ptStar>3) continue;
 
            //
            //insert constituent-level analysis code here
            //
            hPID->Fill((genDau_pid->at(j)).at(k));
            hPIDvsMult->Fill(nmult, (genDau_pid->at(j)).at(k));

            if(fabs((genDau_pid->at(j)).at(k)) == 22)
            {
              if( fabs((genDau_mom->at(j)).at(k))==21 || fabs((genDau_mom->at(j)).at(k))<=6 )
              {
                hPtFrag->Fill(ptStar);
                hPtFragvsMult->Fill(nmult,ptStar);
                hPtBvsPtJFrag->Fill(ptStar, pt);
                hEtaBvsEtaJFrag->Fill(etaStar,eta);
              }

//              hPtBvsPtJFrag->Fill(ptStar, pt);

              hPt->Fill(ptStar);
              hPtvsMult->Fill(nmult,ptStar);
              hMomPID->Fill((genDau_mom->at(j)).at(k));
              hMomPIDvsMult->Fill(nmult,(genDau_mom->at(j)).at(k));
            }
          } //end of constituent loop
        }
      }  //end of jet loop
    }  //end of event loop
    inputFile->Close();
  } //end of file loop
  
  TFile * outputFile = TFile::Open(outputFileName.data(),"recreate");
  hJetPvsMult->Write();
  hPt->Write();
  hPtFrag->Write();
  hPtBvsPtJFrag->Write();
  hEtaBvsEtaJFrag->Write();
  hPtvsMult->Write();
  hPtFragvsMult->Write();
  hPID->Write();
  hPIDvsMult->Write();
  hMomPID->Write();
  hMomPIDvsMult->Write();
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
