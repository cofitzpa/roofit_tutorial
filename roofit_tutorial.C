/* RooFit tutorial
 * 
 * Highlights some of the basic features of RooFit by making 
 * fits to the data provided in dataset.root, increasing in complexity
 *
 * This tutorial written by Conor Fitzpatrick
 * conor.fitzpatrick@cern.ch
 * with thanks to Wouter Verkerke and David Kirkby 
 * If you spot any bugs or have any problems working through this, let me know... 
 * 
 */

#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"

using namespace RooFit;

void roofit_tutorial(){

	/************* SOME BASICS FIRST ***************
	 * First things first, make sure you have access to ROOT. 
	 * To do this on lxplus, type "SetupProject Gaudi ROOT"
	 *
	 * Next, make sure you've got a folder containing this 
	 * macro and the dataset.root ntuple. 
	 *
	 * Open dataset.root in a TBrowser and have a look at 
	 * the ntuple. In it you'll see two columns labelled mass, time. 
	 * This is the same as any other ntuple, there's nothing RooFit specific about it. 
	 *
	 * The first few bits of code in this tutorial will load the dataset from 
	 * ntuple and plot them in the manner of RooFit. This is to get you started. 
	 * Once you know what is going on in the first few lines, 
	 * run this macro in the same directory as the one containing dataset.root:
	 * "root -x roofit_tutorial.C"
	 * If you see the plots produced, you're ready to read on...
	 */

	//First we open the input dataset.root file:
	TFile* inputFile = TFile::Open("dataset.root");
	//Now we attach to the TTree containing the data:
	TTree* inputTree = (TTree*)inputFile->Get("modelData");

	//Define the observables and ranges over which the PDFs will be made 
	//based on the two columns in the ntuple:
	RooRealVar obs_mass("mass","mass",5300,5400,"MeV/c^{2}");
	RooRealVar obs_time("time","time",0.0,12.0,"ps");

	//And create a RooDataSet object containing these observables:
	RooDataSet *data = new RooDataSet("data","data",inputTree,RooArgList(obs_mass,obs_time));

	//Next, let's plot the RooDataSet in each dimension to see what 
	//the distributions look like:

	//Create a TCanvas
	TCanvas *c = new TCanvas("data","data",1024,512);
	//Chop it in two to show both dimensions:
	c->Divide(2);

	//First we'll plot the mass on the left hand side
	c->cd(1);
	//Adjust the margins as the root default is terrible
	gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15);
	//Now the important bit: A RooPlot is made in the mass dimension:
	RooPlot *massplot = obs_mass.frame();

	//We plot the data on this RooPlot:
	data->plotOn(massplot);
	massplot->GetYaxis()->SetTitleOffset(1.6);
	massplot->Draw();

	//Now we plot the time distribution on the RHS: 
	c->cd(2);
	gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15);
	//Logarithmic axis as this looks better:
	gPad->SetLogy();
	RooPlot *timeplot = obs_time.frame();
	data->plotOn(timeplot);
	timeplot->GetYaxis()->SetTitleOffset(1.6);
	timeplot->Draw();
	c->Draw();
	c->Print("dataplot.pdf");

	/***********************TUTORIAL STARTS HERE************************
	 * Your task is to try and fit to this dataset to extract the signal 
	 * yield and lifetime.
	 * We'll start with the signal yield, fitting to the mass distribution.
	 * Then we'll try to extract the signal lifetime from the time distribution
	 * Finally, we'll simultaneously fit to both time and mass to see if it 
	 * improves the errors. 
	 *
	 * The tutorial is written in such a way that you should uncomment 
	 * lines one at a time and fill them in as you go along. This way 
	 * the code can be re-run at each step to see what happens. 
	 */


	/**********************THE MASS FIT********************************
	 * Some hints: The signal mass is a simple gaussian, you can use a 
	 * RooGaussian PDF for the signal model. The background is an O(1) 
	 * polynomial with a negative gradient. Polynomials are 
	 * surprisingly hard to fit, so to help you out I'll tell you that 
	 * the zero order term is 5500 and should be _fixed_, allowing only 
	 * the gradient to float. Limits are important: Is it worth trying 
	 * to fit the first order term on a range >0? Try the range -2.0 -> 0.0
	 */

	//THE SIGNAL MASS PDF
	//To start you out, knowing that the signal PDF is a gaussian 
	//you're going to need two variables to specify it: 
	//A mean and a width. Uncomment these and complete the second line: 

//	RooRealVar signal_mass_mean("signal_mass_mean","signal_mass_mean",5300,5400,"MeV/c^{2}");
//	RooRealVar signal_mass_sigma(

	//Now you'll need to specify the RooGaussian PDF. 
	//Have a look at http://root.cern.ch/root/html/RooGaussian.html 
	//to see what it expects
	//Complete this line, using the obs_mass observable and the 
	//mean, width you just finished off: 

//	RooAbsPdf* signal_mass = new RooGaussian(

	//THE BACKGROUND MASS PDF
	//It's now up to you to create the O(1) poly to fit the background with.
	//First you'll need to specify the coefficients for both O(0) and O(1)
	//Remember that the O(0) term should be _fixed_ and not floated in the fit.

//	RooRealVar bkg_mass_p0(
//	RooRealVar bkg_mass_p1(

	//Now finish this line off with a RooPolynomial (http://root.cern.ch/root/html/RooPolynomial.html)

//	RooAbsPdf* bkg_mass = new RooPolynomial(


	//THE MASS FIT
	// Now we have a PDF describing the signal and another describing 
	// the background in the mass observable. We can combine these as 
	// a RooAddPdf with an additional 2 parameters describing the yields.
	// When we fit, assuming the PDFs are a good description of the data 
	// we will get back the signal mean, the width and the yield. 


	//We need two RooRealVars to specify the signal and background yields. 
	// As we don't know the yields a priori, we'll set the ranges for each 
	// to be from 0 to the number of events in the ntuple. Uncomment these:

//	RooRealVar Nsig("Nsig","Nsig",0.0,data->numEntries());
//	RooRealVar Nbkg("Nbkg","Nbkg",0.0,data->numEntries());

	//Now we create a RooAddPdf, telling it to use these yield terms 
	// and the signal and background pdfs in the same order: 
//	RooAbsPdf *mass_pdf = new RooAddPdf("mass_pdf","mass_pdf",RooArgList(*signal_mass,*bkg_mass),RooArgList(Nsig,Nbkg));

	//Finally, the extended likelihood fit: We fit to the mass_pdf and save the fitresult: 
//	RooFitResult *mass_result = mass_pdf->fitTo(*data,Extended(),Save());
	
	//We print the mass fit result to the screen: 
//	mass_result->Print();

	//We update the mass plot from earlier to show the fit to the data. 
	//Uncommenting the next few lines will plot the PDF and its components:
//	c->cd(1);
	//The combined fit in solid blue:	
//	mass_pdf->plotOn(massplot);
	//The signal component only in dashed blue:
//	mass_pdf->plotOn(massplot,Components(*signal_mass),LineColor(kBlue),LineStyle(kDashed));
	//The background component only in dashed red: 
//	mass_pdf->plotOn(massplot,Components(*bkg_mass),LineColor(kRed),LineStyle(kDashed));
//	massplot->Draw();
//	c->Update();
	// Save a PDF of current plot:
//	c->Print("masspdfplot.pdf");

	//Once you get to here, make sure the code runs and that the plot 
	//gets updated with the mass fit. Does it look reasonable? 
	//Look at the printed fit result. How many signal events do you think I generated? 
	//It's a nice round number....
	//If all has worked, proceed to the lifetime fit. 
	     

	//******************** THE LIFETIME FIT ***************************
	/* Now we'll try to fit the lifetime of the signal distribution
	 * RooFit comes with a few pre-built "B-decay specific" PDFs that 
	 * include CP violation, resolution modelling, mistag rates, etc. 
	 * These are a little too overkill for a basic tutorial, so we'll 
	 * concern ourselves with a simple exponential decay model here. 
	 * 
	 * Part of the power of RooFit is the way in which PDFs can be 
	 * modified using formulae. The RooExponential PDF is of the form
	 * f = e^{c*x} but we want f = e^{-x/c}. We can do this with 
	 * a RooFormulaVar
	 */

	// THE SIGNAL TIME PDF
	// First, we specify the lifetime parameter as a RooRealVar 
	// with a range as usual:
//	RooRealVar signal_lifetime("signal_lifetime","signal_lifetime",0.0,5.0,"ps");

	//Now we need to turn this into a derived variable of the 
	// type the PDF expects, in this case c = -1/tau. 
	// See http://root.cern.ch/root/html/RooFormulaVar.html and 
	// complete the line below: 
//	RooFormulaVar signal_exponent(

	//Finally we give the PDF this term as the free parameter in the exponential PDF
	//Complete this line:
//	RooAbsPdf* signal_time = new RooExponential(

	//THE BACKGROUND TIME PDF
	//The background in this case is "prompt"- It has a much smaller lifetime 
	//than the signal but can be considered as roughly exponential. 
	//We make an identical PDF to the signal one. Uncomment and complete:
	
//	RooRealVar bkg_lifetime("bkg_lifetime","bkg_lifetime",0.0,3.0,"ps");
//	RooFormulaVar bkg_exponent(
//	RooAbsPdf* bkg_time = new RooExponential(

	//THE LIFETIME FIT
	// Now we can do the exact same thing as we did with the mass PDF, 
	// and fit for the lifetimes

	// We don't actually need new RooRealVars to specify the yields, 
	// we can recycle the old ones from the mass fit
	// Uncomment and complete the time PDF and the fit, being sure to save the fitresult: 
//	RooAbsPdf *time_pdf = new RooAddPdf(
//	RooFitResult *time_result = time_pdf->fitTo(
//	time_result->Print();

	//Now we update the time plot as we did for the mass:
//	c->cd(2);
	//The combined fit in solid blue:
//	time_pdf->plotOn(timeplot);
	//The signal component only in dashed blue:
//	time_pdf->plotOn(timeplot,Components(*signal_time),LineColor(kBlue),LineStyle(kDashed));
	//The background component only in dashed red: 
//	time_pdf->plotOn(timeplot,Components(*bkg_time),LineColor(kRed),LineStyle(kDashed));
//	timeplot->Draw();
//	c->Update();
//	c->Print("timepdfplot.pdf");

	//Again, now that you've gotten this far, make sure the code compiles and that 
	//the fit result looks good. 
	//If you know your particles the mass + lifetime result should look familiar :) 
	//You can now go on to the 2D fit. 

	/******************* THE FULL 2D PDF *****************
	 * Finally we combine the 2 PDFs to make a 2D RooProductPDF 
	 * in the 2 observables
	 * By fitting to both dimensions simultaneously we make more information available 
	 * to the fit which leads to a decrease in the size of the errors on the fit parameters
	 */

	//We make a product PDF of the two signal PDFs first. This 2D PDF completely 
	//specifies the signal in time and mass:
//	RooAbsPdf* signal = new RooProdPdf(

	//Likewise for the background:
//	RooAbsPdf* bkg = new RooProdPdf(


	// And now we make an AddPdf of the signal + background as before, 
	// but this time in both dimensions using the 2D PDFs:
//	RooAbsPdf* model = new RooAddPdf(

	//Now we fit to the full model:
//	RooFitResult *model_result = model->fitTo(*data,Extended(),Save());

	//Because we saved the fit results, we can print these sequentially and compare: 
//	cout << "------------------ FIT RESULT FOR MASS ONLY --------------" << endl;
//	mass_result->Print();

//	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
//	time_result->Print();

//	cout << "----------------- FIT RESULT FOR THE 2D MODEL ------------" << endl;
//	model_result->Print();
 	
	//If you've gotten this far, make sure the code runs to completion, and read 
        //through the output. Compare the printed fit results at the very end. You 
        //should see a decrease in the error on all parameters in the 2D model. 
        //If you do and still have time, read on to see a toy study in action....
 

	//****************************TOY STUDIES*********************************
	/* Lastly, we do a short toy study based on the result of our fit: What is 
	 * the _expected_ sensitivity of our model to the lifetime assuming the 
	 * result we measured? In my opinion this is the coolest feature of RooFit: 
	 * The models you build can be fit to data, can generate toy data or can 
	 * run big toy studies only with a couple of extra lines of code..
	 * This last block requires no modification to run- just remove the 
	 * slash-star and star-slash to uncomment it, then re-run the code.
	 */

/*
	RooMCStudy *toyStudy = new RooMCStudy(*model,RooArgList(obs_mass,obs_time),Extended(),FitOptions(Save(kTRUE)),Silence());
	UInt_t Ntoys = 50;
	toyStudy->generateAndFit(Ntoys);
	// That's it! Only 2 lines to automatically generate and fit Ntoys data 
	// samples with the parameters of the model fit result. 
	// Do you understand what the above lines are doing? 
	
	//Now that we've run our toys, we want to see what our expected 
	//sensitivity to the lifetime is. We'll plot the signal lifetime, 
	//its error and pull for the toys:  
	
	TCanvas *d = new TCanvas("toys","toys",1024,340);
	d->Divide(3);
	d->cd(1);
	RooPlot *value = toyStudy->plotParam(signal_lifetime,Bins(10));
	value->Draw();
	d->cd(2);
	RooPlot *error = toyStudy->plotError(signal_lifetime,Bins(10));
	error->Draw();
	d->cd(3);
	RooPlot *pull = toyStudy->plotPull(signal_lifetime,Bins(10),FitGauss(kTRUE));
	pull->Draw();
	d->Update();
	d->Print("toystudy1.pdf");

	//We can also extract the fit result for each individual fit, 
	//and look at correlations between fitted values:
	
	TCanvas *e = new TCanvas("toys2","toys2",1280,384);
	e->Divide(3);
	e->cd(1);
	gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15) ;
	
	// A 2D scatter showing how the signal lifetime and background lifetime fluctuate
	TH1* sig_vs_bkg_lifetime = toyStudy->fitParDataSet().createHistogram("sig_vs_bkg_lifetime",signal_lifetime,YVar(bkg_lifetime)) ;
	sig_vs_bkg_lifetime->SetTitle("Signal vs Background Lifetime");
	sig_vs_bkg_lifetime->Draw("box");

	e->cd(2);
	gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15) ;
	
	// Similarly for the yields:
	TH1* sig_vs_bkg_yield = toyStudy->fitParDataSet().createHistogram("sig_vs_bkg_yield",Nsig,YVar(Nbkg));
	sig_vs_bkg_yield->SetTitle("Signal vs Background Yield");
	sig_vs_bkg_yield->Draw("box");

	e->cd(3);
	gPad->SetLeftMargin(0.20); gPad->SetBottomMargin(0.15);

	// We extract the correlation matrix of each toy and average them to get the average correlation matrix across 50 toys:
	TH2* corrmatrix = toyStudy->fitResult(0)->correlationHist("Parameter Correlation matrix");
	for(UInt_t i=1; i<Ntoys; i++){
		          TString name = "c";
			            name += i;
				      corrmatrix->Add(toyStudy->fitResult(i)->correlationHist(name));
	}
	corrmatrix->Scale(1.0/((double)Ntoys));
	corrmatrix->SetMinimum(-1.0);

	gStyle->SetOptStat(0);
	corrmatrix->Draw("colz");
	e->Draw();
	e->Print("toystudy2.pdf");
*/

	//Here ends the tutorial. If you want to go further with RooFit, the tutorials 
	//bundled in ROOT are reasonably good. On lxplus you can find these in the 
	// $ROOTSYS/tutorials/roofit/ directory  
	// There is a lot more that RooFit can do. I hope this tutorial gave you a feel for the basics!

}
