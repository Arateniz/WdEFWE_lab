#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <vector>
#include <cmath>

struct TLV {
	TLorentzVector v1;
	TLorentzVector v2;
	TLorentzVector v3;
	double w;
};

#ifdef __ROOTCLING__
#pragma link C++ class ParticlePair+;
#pragma link C++ class std::vector<ParticlePair>+;
#endif

template <typename T>
std::pair<int, int> getTwo(
	const std::vector<T> &vec
){
/* vim-marker---> */
	if (vec.size() < 2) {
		throw std::invalid_argument("Vector length < 2.");
	}

	T max1 = -std::numeric_limits<T>::infinity();
	T max2 = -std::numeric_limits<T>::infinity();

	int leading    = 0;
	int subleading = 0;
	int i = 0;

	for (double val : vec) {
		if (val > max1) {
			max2 = max1;
			max1 = val;
			leading = i;
		} else if (val > max2) {
			max2 = val;
			subleading = i;
		}
		++i;
	}
/* <---vim-marker */
	return {leading, subleading};
}

void selection(
	std::initializer_list<TTree *> trees,
	std::vector<TLV> *out,
	double lumi,
	TH1F ** cfHist = nullptr,
	const std::string &histName = "hist"
){
	/* vim-marker---> */
	int nTrees = trees.size();
	int nTree  = 0;

	if (cfHist && *cfHist == nullptr) {
		*cfHist = new TH1F(
			histName.c_str(),
			"Cut Flow;Selection Step;Events",
			8, 0.5, 8.5
		);

		(*cfHist)->GetXaxis()->SetBinLabel(1, "Total");
		(*cfHist)->GetXaxis()->SetBinLabel(2, "TightID");
		(*cfHist)->GetXaxis()->SetBinLabel(3, "pT");
		(*cfHist)->GetXaxis()->SetBinLabel(4, "Eta");
		(*cfHist)->GetXaxis()->SetBinLabel(5, "Crack veto");
		(*cfHist)->GetXaxis()->SetBinLabel(6, "Iso 1");
		(*cfHist)->GetXaxis()->SetBinLabel(7, "Iso 2");
		(*cfHist)->GetXaxis()->SetBinLabel(8, "All Cuts");
	}

	for (auto tree : trees) {
		nTree++;
		std::pair<int, int> p = {0, 0};

		std::vector<int> *isTight = nullptr;
		std::vector<double> *pt30 = nullptr;
		std::vector<double> *et20 = nullptr;

		std::vector<double> *pt  = nullptr;
		std::vector<double> *eta = nullptr;
		std::vector<double> *phi = nullptr;
		std::vector<double> *E   = nullptr;

		float Xsec  = 0.0;
		float SWght = 0.0;
		float MWght = 0.0;

		double weight = 0.0;

		tree->SetBranchAddress("photon_isTightID", &isTight);
		tree->SetBranchAddress("photon_ptcone30",  &pt30);
		tree->SetBranchAddress("photon_etcone20",  &et20);

		tree->SetBranchAddress("photon_pt",  &pt);
		tree->SetBranchAddress("photon_eta", &eta);
		tree->SetBranchAddress("photon_phi", &phi);
		tree->SetBranchAddress("photon_E",   &E);

		tree->SetBranchAddress("XSection",   &Xsec);
		tree->SetBranchAddress("SumWeights", &SWght);
		tree->SetBranchAddress("mcWeight",   &MWght);

		int last_percent = -1;
		Long64_t nEntries = tree->GetEntries();
		for (Long64_t i = 0; i < nEntries; ++i) {
			int percent = static_cast<int>(100.0 * i / nEntries);
			if (percent != last_percent) {
				std::cout << "\rSelection "
					<< nTree
					<< " of "
					<< nTrees
					<< ", progress: "
					<< percent
					<< "%"
					<< std::flush;

				last_percent = percent;
			}

			tree->GetEntry(i);
			p = getTwo<double>(*pt);

			MWght  = MWght==0?1:MWght;
			weight = lumi * Xsec / SWght * MWght;

			if (cfHist && *cfHist) (*cfHist)->Fill(1, weight);

			bool cut1 = (
				(isTight->at(p.first))
					&&
				(isTight->at(p.second))
			);
			if (!cut1) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(2, weight);

			bool cut2 = (
				(pt->at(p.first) > 35e3)
					&&
				(pt->at(p.second) > 25e3)
			);
			if (!cut2) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(3, weight);

			bool cut3 = (
				(std::fabs(eta->at(p.first)) < 2.37)
					&&
				(std::fabs(eta->at(p.second)) < 2.37)
			);
			if (!cut3) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(4, weight);

			bool cut4 = (
				(std::fabs(eta->at(p.first)) > 1.52)
					||
				(std::fabs(eta->at(p.first)) < 1.37)
			);
			if (!cut4) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(5, weight);

			bool cut5 = (
				(std::fabs(eta->at(p.second)) > 1.52)
					||
				(std::fabs(eta->at(p.second)) < 1.37)
			);
			if (!cut5) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(6, weight);

			bool cut6 = (
				(pt30->at(p.first)/pt->at(p.first) < 0.065)
					&&
				(et20->at(p.first)/pt->at(p.first) < 0.065)
			);
			if (!cut6) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(7, weight);

			bool cut7 = (
				(pt30->at(p.second)/pt->at(p.second) < 0.065)
					&&
				(et20->at(p.second)/pt->at(p.second) < 0.065)
			);
			if (!cut7) continue;
			if (cfHist && *cfHist) (*cfHist)->Fill(8, weight);

			TLorentzVector vec1, vec2;

			vec1.SetPtEtaPhiE(
				pt->at(p.first)/1e3,
				eta->at(p.first),
				phi->at(p.first),
				E->at(p.first)/1e3
			);

			vec2.SetPtEtaPhiE(
				pt->at(p.second)/1e3,
				eta->at(p.second),
				phi->at(p.second),
				E->at(p.second)/1e3
			);

			out->push_back((TLV){vec1, vec2, vec1 + vec2, weight});
		}

		std::cout << "\r\rSelection "
			<< nTree
			<< " of "
			<< nTrees
			<< ", progress: 100%"
			<< std::endl;
	}
	/* <---vim-marker */
}

template<typename T>
void genericHist(
	std::initializer_list<TTree *> trees,
	const std::string &dir,
	const std::string &title,
	const std::string &var_name,
	const std::string &unit,
	const int bins,
	const double lbound,
	const double rbound,
	const bool splitleading = false,
	const bool leading = true
){
/* vim-marker---> */
	TH1F *hist = new TH1F(
		var_name.c_str(),
		(title + ";" + var_name + " ["+unit+"];Counts").c_str(),
		bins, lbound, rbound
	);

	if (splitleading) {
		if (leading)
			hist->GetXaxis()->SetTitle((var_name + " leading" + " ["+unit+"]").c_str());
		else
			hist->GetXaxis()->SetTitle((var_name + " subleading" + " ["+unit+"]").c_str());
	}


	double scale = 1.0;
	for (TTree* tree : trees) {
		tree->ResetBranchAddresses();
		if constexpr (std::is_same<T, float>::value) {
			std::vector<T>* val_vector = nullptr;
			tree->SetBranchAddress(var_name.c_str(), &val_vector);
			if (unit == "GeV") scale = 1e-3;

			Long64_t nEntries = tree->GetEntries();
			for (Long64_t i = 0; i < nEntries; ++i) {
				tree->GetEntry(i);
				if (!val_vector) continue;

				if (splitleading) {
					if (val_vector->size() < 2) continue;
					auto p = getTwo<T>(*val_vector);
					hist->Fill(leading ? (*val_vector)[p.first]*scale:(*val_vector)[p.second]*scale);
				} else {
					for (const auto& val : *val_vector) {
						hist->Fill(val*scale);
					}
				}
			}
		} else {
			T val_scalar;
			tree->SetBranchAddress(var_name.c_str(), &val_scalar);

			Long64_t nEntries = tree->GetEntries();
			for (Long64_t i = 0; i < nEntries; ++i) {
				tree->GetEntry(i);
				hist->Fill(val_scalar);
			}
		}
	}

	TCanvas *c = new TCanvas(
		var_name.c_str(),
		title.c_str(),
		600, 500
	);

	c->cd();
	gPad->SetLeftMargin(0.15);
	hist->Draw("hist");

	c->Update();
	TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
	if (stats) {
		stats->SetX1NDC(0.70);
		stats->SetX2NDC(0.90);
		stats->SetY1NDC(0.75);
		stats->SetY2NDC(0.90);
		stats->SetTextSize(0.03);
	}

	c->Modified();
	c->Update();


	std::string filename = dir + var_name + ".pdf";

	if (splitleading) {
		if (leading) filename = dir + var_name + "-leading.pdf";
		else filename = dir + var_name + "-subleading.pdf";
	}

	c->SaveAs(filename.c_str());

	delete c;
	delete hist;
/* <---vim-marker */
}

void drawHist(
	std::vector<TLV> *vec,
	const std::string &dir,
	const std::string &lep_name,
	const std::string &var_name,
	const std::string &unit,
	const int bins,
	const double lbound,
	const double rbound,
	std::function<double(const TLV&)> eval,
	const std::vector<TF1*> &fitFunctions = {},
	bool coverup = false,
	double lcvr = std::numeric_limits<double>::quiet_NaN(),
	double rcvr = std::numeric_limits<double>::quiet_NaN()
){
	/* vim-marker---> */
	if (coverup && (std::isnan(lcvr) || std::isnan(rcvr))) {
		std::cout << "Invalid cover-up bounds!" << std::endl;
		return;
	}

	TH1F *hist = new TH1F(
		var_name.c_str(),
		(lep_name + ";" + var_name + " [" + unit + "];Counts").c_str(),
		bins, lbound, rbound
	);

	TCanvas *c = new TCanvas(
		var_name.c_str(),
		lep_name.c_str(),
		600, 500
	);

	c->cd();
	gPad->SetLeftMargin(0.15);

	for (Long64_t i = 0; i < vec->size(); ++i) {
		double val = eval(vec->at(i));
		if (coverup && val > lcvr && val < rcvr) continue;
		hist->Fill(val, vec->at(i).w);
	}

	if (!std::isnan(lcvr) && !std::isnan(rcvr)) {
		int bin_min = hist->GetXaxis()->FindBin(120);
		int bin_max = hist->GetXaxis()->FindBin(130);
		double sum  = hist->Integral(bin_min, bin_max);

		std::cout << "Integral of histogram '" << hist->GetName()
		<< "' in range [" << lcvr << ", " << rcvr << "] = "
		<< sum << std::endl;
	}

	TLegend *legend = new TLegend(0.7, 0.75, 0.9, 0.67);

	hist->Draw("hist");

	bool validFit = false;
	for (TF1* f : fitFunctions) {
		if (f != nullptr) {
			f->SetNpx(1000);
			hist->Fit(f, "R");
			f->Draw("same");
			legend->AddEntry(f, f->GetName());
			validFit = true;
		}
	}

	if (validFit)
		legend->Draw();

	c->Update();

	if (validFit) {
		for (TF1* f : fitFunctions) {
			if (f != nullptr) {
				double integral = f->Integral(lcvr, rcvr);
				std::cout << "Integral of fit function '" << f->GetName()
					<< "' in range [" << lcvr << ", " << rcvr << "] = "
					<< integral << std::endl;
			}
		}
	}

	TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
	if (stats) {
		stats->SetX1NDC(0.70);
		stats->SetX2NDC(0.90);
		stats->SetY1NDC(0.75);
		stats->SetY2NDC(0.90);
		stats->SetTextSize(0.03);
	}

	c->Modified();
	c->Update();

	c->SaveAs((dir + var_name + ".pdf").c_str());

	delete c;
	delete hist;
	delete legend;
	/* <---vim-marker */
}

void drawCutFlow(
	TH1F **hist,
	const std::string &dir,
	const std::string &sample_name,
	const double scale
){
	/* vim-marker---> */
	if (!hist || !*hist)
		return;

	TCanvas *c = new TCanvas(
		"cutflow_canvas", "Cut Flow",
		600, 500
	);

	(*hist)->SetStats(0);
	(*hist)->SetFillColor(kAzure + 7);
	(*hist)->LabelsOption("h", "X");
	(*hist)->Draw("hist");

	for (int i = 1; i <= (*hist)->GetNbinsX(); ++i) {
		double val = (*hist)->GetBinContent(i);

		if (val == 0)
			continue;

		double x = (*hist)->GetXaxis()->GetBinCenter(i);
		double y = val + scale * (*hist)->GetMaximum();

		TLatex *label = new TLatex(x, y, Form("%.0f", val));
		label->SetTextAlign(22);
		label->SetTextSize(0.03);
		label->Draw();
	}

	c->SaveAs((dir + "Higgs-" + sample_name + "-CutFlow.pdf").c_str());

	delete c;
	/* <---vim-marker */
}

void drawCompare(
	std::initializer_list<std::vector<TLV>*> vectors,
	const std::string &dir,
	const std::string &lep_name,
	const std::string &var_name,
	const std::string &unit,
	const int bins,
	const double lbound,
	const double rbound,
	std::function<double(const TLV&)> eval,
	const std::vector<TF1*>& fitFunctions = {}
){
/* vim-marker---> */
	auto it = vectors.begin();
	std::vector<TLV>* exVec = *it++;
	std::vector<TLV>* mcVec = *it++;

	TH1F* h_data = new TH1F("h_data", "", bins, lbound, rbound);
	TH1F* h_mc   = new TH1F("h_mc", "", bins, lbound, rbound);

	for (const TLV &tlv : *exVec)
		h_data->Fill(eval(tlv), tlv.w);

	for (const TLV &tlv : *mcVec)
		h_mc->Fill(eval(tlv), tlv.w);

	h_mc->SetFillColor(kAzure + 7);
	h_mc->SetLineColor(kBlack);

	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.9);
	h_data->SetLineColor(kBlack);

	TCanvas* c = new TCanvas("c", "canvas", 800, 600);
	c->Divide(1, 2, 0.0, 0.0);

	c->cd(1);
	gPad->SetPad(0, 0.3, 0.95, 0.95);
	gPad->SetBottomMargin(0);
	gStyle->SetOptStat(0);

	TH1F* hFrame = new TH1F("hFrame", "", bins, lbound, rbound);
	hFrame->SetTitle((lep_name + "; " + var_name + " [" + unit + "]; Counts").c_str());
	hFrame->SetMinimum(0);
	hFrame->SetMaximum(std::max(h_data->GetMaximum(), h_mc->GetMaximum()) * 1.4);
	hFrame->GetYaxis()->SetTitleOffset(1.2);
	hFrame->Draw();

	h_mc->Draw("hist same");
	h_data->Draw("E1 same");

	for (TF1* func : fitFunctions) {
		func->SetLineWidth(2);
		func->Draw("same");
	}

	TLegend* legend = new TLegend(0.7, 0.75, 1, 1);
	legend->AddEntry(h_data, "Data", "lep");
	legend->AddEntry(h_mc, "sig", "f");

	for (TF1* func : fitFunctions)
		legend->AddEntry(func, func->GetTitle(), "l");

	legend->Draw();

	c->cd(2);
	gPad->SetPad(0, 0.0, 0.95, 0.3);
	gPad->SetTopMargin(0.0);
	gPad->SetBottomMargin(0.3);

	TH1F* h_diff = (TH1F*) h_data->Clone("h_diff");
	h_diff->Reset();

	if (fitFunctions.size() >= 1) {
		TF1* bkg = fitFunctions[1];
		for (int i = 1; i <= h_data->GetNbinsX(); ++i) {
			double x = h_data->GetBinCenter(i);
			double dataVal = h_data->GetBinContent(i);
			double bkgVal  = bkg->Eval(x);
			h_diff->SetBinContent(i, dataVal - bkgVal);
			h_diff->SetBinError(i, h_data->GetBinError(i));
		}
	}

	h_diff->SetTitle(("; " + var_name + " [" + unit + "]; Data - bkg").c_str());
	h_diff->GetYaxis()->SetTitleSize(0.08);
	h_diff->GetYaxis()->SetTitleOffset(0.5);
	h_diff->GetYaxis()->SetLabelSize(0.08);
	h_diff->GetXaxis()->SetLabelSize(0.08);
	h_diff->GetXaxis()->SetTitleSize(0.1);
	h_diff->GetXaxis()->SetTitleOffset(1.0);
	h_diff->SetMarkerStyle(20);
	h_diff->SetMarkerSize(0.9);
	h_diff->Draw("E1");

	int bin_min = h_diff->GetXaxis()->FindBin(120);
	int bin_max = h_diff->GetXaxis()->FindBin(130);
	double sum  = h_diff->Integral(bin_min, bin_max);

	std::cout << "Integral of histogram '" << h_diff->GetName()
	<< "' in range [" << lbound << ", " << rbound << "] = "
	<< sum << std::endl;

	if (fitFunctions.size() >= 2) {
		TF1* sig = fitFunctions[0];
		TF1* bkg = fitFunctions[1];

		std::string name = "f_diff_sig_bkg";
		TF1* diff = new TF1(
			name.c_str(),
			[sig, bkg](double* x, double*) {
				return sig->Eval(x[0]) - bkg->Eval(x[0]);
			},
			lbound, rbound, 0
		);

		diff->SetLineColor(kRed);
		diff->SetLineStyle(1);
		diff->SetLineWidth(2);
		diff->SetTitle("Signal - Bkg");
		diff->Draw("same");
	}

	TLine* zero = new TLine(lbound, 0, rbound, 0);
	zero->SetLineStyle(2);
	zero->SetLineColor(kAzure);
	zero->Draw("same");

	c->SaveAs((dir + var_name + ".pdf").c_str());

	delete c;
/* <---vim-marker */
}

Double_t poly3(Double_t *x, Double_t *par)
{
	return par[0] + par[1]*x[0] + par[2]*std::pow(x[0], 2) + par[3]*std::pow(x[0], 3);
}

void prog()
{
	TFile *monte_tt = TFile::Open("data/mc_341081.ttH125_gamgam.GamGam.root");
	TFile *monte_gg = TFile::Open("data/mc_343981.ggH125_gamgam.GamGam.root");
	TFile *monte_vb = TFile::Open("data/mc_345041.VBFH125_gamgam.GamGam.root");
	TFile *monte_wp = TFile::Open("data/mc_345318.WpH125J_Wincl_gamgam.GamGam.root");
	TFile *monte_zh = TFile::Open("data/mc_345319.ZH125J_Zincl_gamgam.GamGam.root");

	TTree *tree_tt = (TTree*)monte_tt->Get("mini;1");
	TTree *tree_gg = (TTree*)monte_gg->Get("mini;1");
	TTree *tree_vb = (TTree*)monte_vb->Get("mini;1");
	TTree *tree_wp = (TTree*)monte_wp->Get("mini;1");
	TTree *tree_zh = (TTree*)monte_zh->Get("mini;1");

	TFile *data_A = TFile::Open("data/data_A.GamGam.root");
	TFile *data_B = TFile::Open("data/data_B.GamGam.root");
	TFile *data_C = TFile::Open("data/data_C.GamGam.root");
	TFile *data_D = TFile::Open("data/data_D.GamGam.root");

	TTree *tree_A = (TTree*)data_A->Get("mini;1");
	TTree *tree_B = (TTree*)data_B->Get("mini;1");
	TTree *tree_C = (TTree*)data_C->Get("mini;1");
	TTree *tree_D = (TTree*)data_D->Get("mini;1");

	genericHist<UInt_t>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "lep_n", "-",
		2, 0, 2
	);

	genericHist<UInt_t>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "photon_n", "-",
		5, 0, 5
	);

	genericHist<UInt_t>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "jet_n", "-",
		13, 0, 10
	);

	genericHist<float>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "photon_pt", "GeV",
		30, 0, 300
	);

	genericHist<float>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "photon_eta", "-",
		30, -3, 3
	);

	genericHist<float>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "photon_phi", "rad",
		30, -15, 15
	);

	genericHist<float>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "photon_pt", "GeV",
		30, 0, 300,
		true, true
	);

	genericHist<float>(
		{tree_gg},
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "photon_pt", "GeV",
		30, 0, 300,
		true, false
	);

	std::vector<TLV> *photons_MC = new std::vector<TLV>();
	std::vector<TLV> *photons_EX = new std::vector<TLV>();

	TH1F *cfMCHist = nullptr;
	TH1F *cfEXHist = nullptr;

	std::cout << "Monte-Carlo:" << std::endl;
	selection(
		{tree_tt, tree_gg, tree_vb, tree_wp, tree_zh},
		photons_MC, 10e3, &cfMCHist, "cfMC"
	);

	std::cout << std::endl << "Experimental data:" << std::endl;
	selection(
		{tree_A, tree_B, tree_C, tree_D},
		photons_EX, 1, &cfEXHist, "cfEX"
	);

	drawCutFlow(&cfMCHist, "pix/", "MC", 0.01);
	drawCutFlow(&cfEXHist, "pix/", "EX", 0.02);

	double bins   = 30;
	double lbound = 100;
	double rbound = 160;

	TF1 *gauss     = new TF1("sig", "gaus", lbound, rbound);
	TF1 *poly3_bkg = new TF1("bkg", poly3, lbound, rbound, 4);
	TF1 *polygauss = new TF1(
		"sig + bkg",
		[poly3_bkg, gauss](double *x, double *p) {
			double xv = x[0];
			double val1 = poly3_bkg->Eval(xv);
			double val2 = gauss->Eval(xv);
			return val1 + val2;
		},
		lbound, rbound, 0
	);

	poly3_bkg->SetLineStyle(2);
	poly3_bkg->SetLineColor(kAzure);
	polygauss->SetLineColor(kRed);

	drawHist(
		photons_MC,
		"pix/Higgs-MC-", "H #rightarrow #gamma#gamma", "Mass", "GeV",
		bins, lbound, rbound,
		[](const TLV &v) {return v.v3.M();},
		{gauss},
		false, 120, 130
	);

	drawHist(
		photons_EX,
		"pix/Higgs-EXcvr-", "H #rightarrow #gamma#gamma", "Mass", "GeV",
		bins, lbound, rbound,
		[](const TLV &v) {return v.v3.M();},
		{poly3_bkg},
		true, 120, 130
	);

	drawCompare(
		{photons_EX, photons_MC},
		"pix/Higgs-Comp-", "H #rightarrow #gamma#gamma", "Mass", "GeV",
		bins, lbound, rbound,
		[](const TLV &v) {return v.v3.M();},
		{polygauss, poly3_bkg}
	);

	delete photons_MC; delete cfMCHist;
	delete photons_EX; delete cfEXHist;

	delete gauss; delete poly3_bkg; delete polygauss;

	gApplication->Terminate();
}
