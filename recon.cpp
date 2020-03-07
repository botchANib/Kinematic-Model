#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TLatex.h"
 
#include <math.h>
#include <string.h>
#include <iostream>

#include "particle.hpp" 
using namespace std;

// This program reconstructs the tau mass and then the b meson mass using the root file "BKstTauMu.root".

void recon() {
	std::cout << "Running..." << std::endl;
    char B_graph = 's', t_graph = 's', n_graph = 's', B_stats = 's', t_stats = 's', n_stats = 's';
	string exi;
	bool checker = true;
	while(checker == true) {
		cout << "Do you want to print the B meson mass?                       (y/N): ";
		cin >> B_graph;
		if(B_graph == 'y' | B_graph == 'Y' | B_graph == 'n' | B_graph == 'N') {
			break;
		} else {
			cout << "You have inputted an incorrect answer, if you would like to exit type 'exit', if you would like to continue type any other letter.\n";
			cin.get();
			getline(cin, exi);
			if(strcmp(exi.c_str(),"exit") == 0) {
				exit(1);
			}
			continue;
		}
	}
	while(checker == true) {
		cout << "Do you want to print the tau lepton mass?                    (y/N): ";
		cin >> t_graph;
		if(t_graph == 'y' | t_graph == 'Y' | t_graph == 'n' | t_graph == 'N') {
			break;
		} else {
			cout << "You have inputted an incorrect answer, if you would like to exit type 'exit', if you would like to continue type any other letter.\n";
			cin.get();
			getline(cin, exi);
			if(strcmp(exi.c_str(),"exit") == 0) {
				exit(1);
			}
			continue;
		}
	}
	while(checker == true) {
		cout << "Do you want to print the tau neutrino 4-momentum components? (y/N): ";
		cin >> n_graph;
		if(n_graph == 'y' | n_graph == 'Y' | n_graph == 'n' | n_graph == 'N') {
			break;
		} else {
			cout << "You have inputted an incorrect answer, if you would like to exit type 'exit', if you would like to continue type any other letter.\n";
			cin.get();
			getline(cin, exi);
			if(strcmp(exi.c_str(),"exit") == 0) {
				exit(1);
			}
			continue;
		}
	}
	while(checker == true) {
		cout << "Do you want to print the B mass statistics?                  (y/N): ";
		cin >> B_stats;
		if(B_stats == 'y' | B_stats == 'Y' | B_stats == 'n' | B_stats == 'N') {
			break;
		} else {
			cout << "You have inputted an incorrect answer, if you would like to exit type 'exit', if you would like to continue type any other letter.\n";
			cin.get();
			getline(cin, exi);
			if(strcmp(exi.c_str(),"exit") == 0) {
				exit(1);
			}
			continue;
		}
	}
	while(checker == true) {
		cout << "Do you want to print the tau mass statistics?                (y/N): ";
		cin >> t_stats;
		if(t_stats == 'y' | t_stats == 'Y' | t_stats == 'n' | t_stats == 'N') {
			break;
		} else {
			cout << "You have inputted an incorrect answer, if you would like to exit type 'exit', if you would like to continue type any other letter.\n";
			cin.get();
			getline(cin, exi);
			if(strcmp(exi.c_str(),"exit") == 0) {
				exit(1);
			}
			continue;
		}
	}
	while(checker == true) {
		cout << "Do you want to print the tau neutrino statistics?            (y/N): ";
		cin >> n_stats;
		if(n_stats == 'y' | n_stats == 'Y' | n_stats == 'n' | n_stats == 'N') {
			break;
		} else {
			cout << "You have inputted an incorrect answer, if you would like to exit type 'exit', if you would like to continue type any other letter.\n";
			cin.get();
			getline(cin, exi);
			if(strcmp(exi.c_str(),"exit") == 0) {
				exit(1);
			}
			continue;
		}
	}
	
	
	
	
	// Inputs to adjust smearing.
	double a, b, c;
	if(B_graph == 'Y' | B_graph == 'y' | B_stats == 'Y' | B_stats == 'y') {
      cout << "For no smearing type 0 to the following:\n";
      cout << "Enter momentum smearing (Recommended 0.005 or 0.01): ";
      cin >> a;
      cout << "Enter x,y position smearing      (Recommended 0.04): ";
      cin >> b;
      cout << "Enter z position smearing         (Recommended 0.2): ";
      cin >> c;
    }
    else {
      a = 0;
      b = 0;
      c = 0;
	}
	
  TFile* tfile = TFile::Open("BKstTauTau.root");
  TTree* ttree = (TTree*) tfile->Get("BKstTauTauTuple/MCDecayTree");

  // Defining 4-momenta of particles using hpp file.
  Particle< Double_t > Kplus_M(    "Kplus",    ttree );
  Particle< Double_t > piminus_M(  "piminus",  ttree );
  Particle< Double_t > tauplus_M(   "tauplus",   ttree );
  Particle< Double_t > piplus_M(   "piplus",   ttree );
  Particle< Double_t > piminus0_M( "piminus0", ttree );
  Particle< Double_t > piminus1_M( "piminus1", ttree );

  // hpp file does not currently work for end vertex positions, so old method being used.
  // Defining B origin vertex components
  double B0_TRUEORIGINVERTEX_X, B0_TRUEORIGINVERTEX_Y, B0_TRUEORIGINVERTEX_Z;
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_X",&B0_TRUEORIGINVERTEX_X);
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_Y",&B0_TRUEORIGINVERTEX_Y);
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_Z",&B0_TRUEORIGINVERTEX_Z);
  // Defining tauminus origin vertex components
  double tauminus_TRUEORIGINVERTEX_X, tauminus_TRUEORIGINVERTEX_Y, tauminus_TRUEORIGINVERTEX_Z;
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_X", &tauminus_TRUEORIGINVERTEX_X);
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_Y", &tauminus_TRUEORIGINVERTEX_Y);
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_Z", &tauminus_TRUEORIGINVERTEX_Z);
  // Defining tauminus end vertex components
  double tauminus_TRUEENDVERTEX_X, tauminus_TRUEENDVERTEX_Y, tauminus_TRUEENDVERTEX_Z;
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_X", &tauminus_TRUEENDVERTEX_X);
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_Y", &tauminus_TRUEENDVERTEX_Y);
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_Z", &tauminus_TRUEENDVERTEX_Z);


  TH1D* hist_Bmass   = new TH1D("B Meson Mass","Reconstructed B Meson Mass",70,0,20000);
  hist_Bmass->SetTitle(";Mass, MeV c^{-2};Number of Entries");
  hist_Bmass->SetFillColor(kRed);
  TH1D* hist_taumass = new TH1D("Tauon Mass","Reconstructed Tau Lepton Meson Mass",80,1450,2100);
  hist_taumass->SetTitle(";Mass, MeV c^{-2};Number of Entries");
  hist_taumass->SetFillColor(kRed);
  TH1D* hist_p9_X    = new TH1D("Tau Neutrino Momentum X","Tau Neutrino Momentum in X",100,0,20000);
  hist_p9_X->SetFillColor(kRed);
  TH1D* hist_p9_Y    = new TH1D("Tau Neutrino Momentum Y","Tau Neutrino Momentum in Y",100,0,20000);
  hist_p9_Y->SetFillColor(kRed);
  TH1D* hist_p9_Z    = new TH1D("Tau Neutrino Momentum Z","Tau Neutrino Momentum in Z",100,0,20000);
  hist_p9_Z->SetFillColor(kRed);
  TH1D* hist_p9_E    = new TH1D("Tau Neutrino Energy","Tau Neutrino Energy",100,0,20000);
  hist_p9_E->SetFillColor(kRed);

  TRandom3* rand = new TRandom3();
  
  for (Long64_t i = 0; i < ttree->GetEntries(); i++ ) {
    ttree->GetEntry(i);
    TLorentzVector Kplus_Mvec    = Kplus_M.getVec();
    TLorentzVector piminus_Mvec  = piminus_M.getVec();
    TLorentzVector tauplus_Mvec   = tauplus_M.getVec();
    TLorentzVector piplus_Mvec   = piplus_M.getVec();
    TLorentzVector piminus0_Mvec = piminus0_M.getVec();
    TLorentzVector piminus1_Mvec = piminus1_M.getVec();
	
    // Particle 2-4 and 6-8 mass in tau or mode change mass!
    double s_Kplus    = rand->Gaus(493.677, 0.013);      // from B J et al (2012) particle listings
    double s_piminus  = rand->Gaus(139.570018, 0.00035); // from C Amsler et al (2008) particle listings
    double s_tauplus  = 1776.86;                       // from Beringer J et al (particle data group) 2012 particle summary
    //double s_muplus = 105.6583755;
    double s_piplus   = rand->Gaus(139.570018, 0.00035);
    double s_piminus0 = rand->Gaus(139.570018, 0.00035);
    double s_piminus1 = rand->Gaus(139.570018, 0.00035);
      
      
      
      
      
      
      
    // Adding errors to particles 2-4
    double Kplus_abs             = sqrt(pow(Kplus_Mvec.X(), 2) + pow(Kplus_Mvec.Y(), 2) + pow(Kplus_Mvec.Z(), 2));
    double Kplus_abs_SMEARED     = rand->Gaus(Kplus_abs, a*Kplus_abs);
    double K_smear_factor        = Kplus_abs_SMEARED/Kplus_abs;
    //double Kplus_X               = rand->Gaus(Kplus_Mvec.X(), K_smear_factor);
    //double Kplus_Y               = rand->Gaus(Kplus_Mvec.Y(), K_smear_factor);
    //double Kplus_Z               = rand->Gaus(Kplus_Mvec.Z(), K_smear_factor);
    double Kplus_X               = (Kplus_Mvec.X()*K_smear_factor);
    double Kplus_Y               = (Kplus_Mvec.Y()*K_smear_factor);
    double Kplus_Z               = (Kplus_Mvec.Z()*K_smear_factor);
    double Kplus_Mvec_E          = sqrt((s_Kplus)*(s_Kplus) + (Kplus_X)*(Kplus_X) + (Kplus_Y)*(Kplus_Y) + (Kplus_Z)*(Kplus_Z));
      
    
    double piminus_abs           = sqrt(pow(piminus_Mvec.X(), 2) + pow(piminus_Mvec.Y(), 2) + pow(piminus_Mvec.Z(), 2));
    double piminus_abs_SMEARED   = rand->Gaus(piminus_abs, a*piminus_abs);
    double piminus_smear_factor  = piminus_abs_SMEARED/piminus_abs;
    double piminus_X             = (piminus_Mvec.X()*piminus_smear_factor);
    double piminus_Y             = (piminus_Mvec.Y()*piminus_smear_factor);
    double piminus_Z             = (piminus_Mvec.Z()*piminus_smear_factor);
    double piminus_Mvec_E        = sqrt((s_piminus)*(s_piminus) + (piminus_X)*(piminus_X) + (piminus_Y)*(piminus_Y) + (piminus_Z)*(piminus_Z));
    
    double tauplus_abs            = sqrt(pow(tauplus_Mvec.X(), 2) + pow(tauplus_Mvec.Y(), 2) + pow(tauplus_Mvec.Z(), 2));
    double tauplus_abs_SMEARED    = rand->Gaus(tauplus_abs, a*tauplus_abs);
    double tauplus_smear_factor   = tauplus_abs_SMEARED/tauplus_abs;
    double tauplus_X              = (tauplus_Mvec.X()*tauplus_smear_factor);
    double tauplus_Y              = (tauplus_Mvec.Y()*tauplus_smear_factor);
    double tauplus_Z              = (tauplus_Mvec.Z()*tauplus_smear_factor);
    double tauplus_Mvec_E         = sqrt((s_tauplus)*(s_tauplus) + (tauplus_X)*(tauplus_X) + (tauplus_Y)*(tauplus_Y) + (tauplus_Z)*(tauplus_Z));
    
    // Adding errors to particles 6-8
    double piplus_abs            = sqrt(piplus_Mvec.X()*piplus_Mvec.X() + piplus_Mvec.Y()*piplus_Mvec.Y() + piplus_Mvec.Z()*piplus_Mvec.Z());
    double piplus_abs_smeared    = rand->Gaus(piplus_abs, a*piplus_abs);
    double piplus_smear_factor   = piplus_abs_smeared/piplus_abs;
    double piplus_X              = (piplus_Mvec.X()*piplus_smear_factor);
    double piplus_Y              = (piplus_Mvec.Y()*piplus_smear_factor);
    double piplus_Z              = (piplus_Mvec.Z()*piplus_smear_factor);
    double piplus_Mvec_E         = sqrt((s_piplus)*(s_piplus) + (piplus_X)*(piplus_X) + (piplus_Y)*(piplus_Y) + (piplus_Z)*(piplus_Z));
      
    double piminus0_abs          = sqrt(piminus0_Mvec.X()*piminus0_Mvec.X() + piminus0_Mvec.Y()*piminus0_Mvec.Y() + piminus0_Mvec.Z()*piminus0_Mvec.Z());
    double piminus0_abs_smeared  = rand->Gaus(piminus0_abs, a*piminus0_abs);
    double piminus0_smear_factor = piminus0_abs_smeared/piminus0_abs;
    double piminus0_X            = (piminus0_Mvec.X()*piminus0_smear_factor);
    double piminus0_Y            = (piminus0_Mvec.Y()*piminus0_smear_factor);
    double piminus0_Z            = (piminus0_Mvec.Z()*piminus0_smear_factor);
    double piminus0_Mvec_E       = sqrt((s_piplus)*(s_piplus) + (piminus0_X)*(piminus0_X) + (piminus0_Y)*(piminus0_Y) + (piminus0_Z)*(piminus0_Z));
      
    double piminus1_abs          = sqrt(piminus1_Mvec.X()*piminus1_Mvec.X() + piminus1_Mvec.Y()*piminus1_Mvec.Y() + piminus1_Mvec.Z()*piminus1_Mvec.Z());
    double piminus1_abs_smeared  = rand->Gaus(piminus1_abs, a*piminus1_abs);
    double piminus1_smear_factor = piminus1_abs_smeared/piminus1_abs;
    double piminus1_X            = (piminus1_Mvec.X()*piminus1_smear_factor);
    double piminus1_Y            = (piminus1_Mvec.Y()*piminus1_smear_factor);
    double piminus1_Z            = (piminus1_Mvec.Z()*piminus1_smear_factor);
    double piminus1_Mvec_E       = sqrt((s_piplus)*(s_piplus) + (piminus1_X)*(piminus1_X) + (piminus1_Y)*(piminus1_Y) + (piminus1_Z)*(piminus1_Z));
      
    // Errors on verticies
    double B0_VERTEX_X = rand->Gaus(B0_TRUEORIGINVERTEX_X, b*sqrt(0.03));
    double B0_VERTEX_Y = rand->Gaus(B0_TRUEORIGINVERTEX_Y, b*sqrt(0.03));
    
    double tauminus_VERTEX_X = rand->Gaus(tauminus_TRUEORIGINVERTEX_X, b);
    double tauminus_VERTEX_Y = rand->Gaus(tauminus_TRUEORIGINVERTEX_Y, b);
    double tauminus_VERTEX_Z = rand->Gaus(tauminus_TRUEORIGINVERTEX_Z, c);

    double tauminus_ENDVERTEX_X = rand->Gaus(tauminus_TRUEENDVERTEX_X, b);
    double tauminus_ENDVERTEX_Y = rand->Gaus(tauminus_TRUEENDVERTEX_Y, b);
    double tauminus_ENDVERTEX_Z = rand->Gaus(tauminus_TRUEENDVERTEX_Z, c);

    // Main Program.
    // Reconstruction calculations.
    double s5unit_X = tauminus_ENDVERTEX_X - tauminus_VERTEX_X;
    double s5unit_Y = tauminus_ENDVERTEX_Y - tauminus_VERTEX_Y;
    double s5unit_Z = tauminus_ENDVERTEX_Z - tauminus_VERTEX_Z;

    double s1unit_X = tauminus_VERTEX_X - B0_VERTEX_X;
    double s1unit_Y = tauminus_VERTEX_Y - B0_VERTEX_Y;
    
    double xi = ((Kplus_X + piminus_X + tauplus_X)*s1unit_Y - (Kplus_Y + piminus_Y + tauplus_Y)*s1unit_X)/((s5unit_Y*s1unit_X) - (s5unit_X*s1unit_Y));
                    
    double p9_X = xi*s5unit_X - (piplus_X + piminus0_X + piminus1_X);
    double p9_Y = xi*s5unit_Y - (piplus_Y + piminus0_Y + piminus1_Y);
    double p9_Z = xi*s5unit_Z - (piplus_Z + piminus0_Z + piminus1_Z);
    double p9_E = sqrt(pow(p9_X, 2) + pow(p9_Y, 2) + pow(p9_Z, 2));
		
    double s5_0 = pow(p9_E + piplus_Mvec_E + piminus0_Mvec_E + piminus1_Mvec_E, 2);
    double s5_1 = pow(p9_X + piplus_X + piminus0_X + piminus1_X, 2) + pow(p9_Y + piplus_Y + piminus0_Y + piminus1_Y, 2) + pow(p9_Z + piplus_Z + piminus0_Z + piminus1_Z, 2);
                    
    double s5 = sqrt(s5_0-s5_1);

    double s4_0 = pow((sqrt(pow(s5, 2) + (pow(s5unit_X, 2) + pow(s5unit_Y, 2) + pow(s5unit_Z, 2))*pow(xi, 2)) + Kplus_Mvec_E + piminus_Mvec_E + tauplus_Mvec_E), 2);

    double s4_1 = pow((s5unit_X*xi + Kplus_X + piminus_X + tauplus_X), 2) + pow((s5unit_Y*xi + Kplus_Y + piminus_Y + tauplus_Y), 2) + pow((s5unit_Z*xi + Kplus_Z + piminus_Z + tauplus_Z), 2);

    double s1 = sqrt((s4_0-s4_1));
    //  cout<<s1<<"\n";
    
    // Histogram being filled with b meson mass values.
    hist_Bmass->Fill(s1);
    hist_taumass->Fill(s5);
    
    hist_p9_X->Fill(p9_X);
    hist_p9_Y->Fill(p9_Y);
    hist_p9_Z->Fill(p9_Z);
    hist_p9_E->Fill(p9_E);
  }
  
  if(B_graph == 'Y' | B_graph == 'y') {
		TCanvas* can1 = new TCanvas("B Meson Mass");
		hist_Bmass->Draw();
	}
	
	if(t_graph == 'Y' | t_graph == 'y') {
		TCanvas* can2 = new TCanvas("Tauon Mass");
		hist_taumass->Draw();
	}
  
  if(n_graph == 'Y' | n_graph == 'y') {
		TCanvas* can3 = new TCanvas("Tau Neutrino Reconstruction");
		can3->Divide(2,2);
		can3->cd(1);
		hist_p9_X->Draw();
		can3->cd(2);
		hist_p9_Y->Draw();
		can3->cd(3);
		hist_p9_Z->Draw();
		can3->cd(4);
		hist_p9_E->Draw();
	}

	
	// Outputs mean and std dev.
  double m_B = hist_Bmass->GetMean();
  double m_tau = hist_taumass->GetMean();
  double m_nX = hist_p9_X->GetMean();
  double m_nY = hist_p9_Y->GetMean();
  double m_nZ = hist_p9_Z->GetMean();
  double m_nE = hist_p9_E->GetMean();
  
  double stdDev_B = hist_Bmass->GetStdDev();
  double stdDev_tau = hist_taumass->GetStdDev();
  double stdDev_nX = hist_p9_X->GetStdDev();
  double stdDev_nY = hist_p9_Y->GetStdDev();
  double stdDev_nZ = hist_p9_Z->GetStdDev();
  double stdDev_nE = hist_p9_E->GetStdDev();

  std::cout << "\n" << "------- Statistics -------" << std::endl;
  
  if(B_stats == 'Y' | B_stats == 'y') {
		std::cout << "\nB mass mean value:         " << m_B << " MeV/c^2" << std::endl;
		std::cout << "B mass standard deviation: " << stdDev_B << " MeV/c^2" << std::endl;
	}
	if(t_stats == 'Y' | t_stats == 'y') {
		std::cout << "\nTauon mass mean value:         " << m_tau << " MeV/c^2" << std::endl;
		std::cout << "Tauon mass standard deviation: " << stdDev_tau << " MeV/c^2" << std::endl;
	}
	if(n_stats == 'Y' | n_stats == 'y') {
		std::cout << "\nTau Neutrino x momentum mean value:         " << m_nX << " MeV/c" << std::endl;
		std::cout << "Tau Neutrino x momentum standard deviation: " << stdDev_nX << " MeV/c" << std::endl;
		std::cout << "\nTau Neutrino y momentum mean value:         " << m_nY << " MeV/c" << std::endl;
		std::cout << "Tau Neutrino y momentum standard deviation: " << stdDev_nY << " MeV/c" << std::endl;
		std::cout << "\nTau Neutrino z momentum mean value:         " << m_nZ << " MeV/c" << std::endl;
		std::cout << "Tau Neutrino z momentum standard deviation: " << stdDev_nZ << " MeV/c" << std::endl;
		std::cout << "\nTau Neutrino energy mean value:         " << m_nE << " MeV/c" << std::endl;
		std::cout << "Tau Neutrino energy standard deviation: " << stdDev_nE << " MeV/c" << std::endl;
	}
	std::cout << "...Finished" << std::endl;
  return;
}

