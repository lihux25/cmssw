#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSimpleRecAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include "RecoLocalCalo/HcalRecAlgos/src/HcalTDCReco.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/rawEnergy.h"

#include <algorithm>
#include <cmath>

#include <TH1.h>
#include <TMinuit.h>

//--- temporary for printouts
// #include<iostream>

constexpr double MaximumFractionalError = 0.002; // 0.2% error allowed from this source

HcalSimpleRecAlgo::HcalSimpleRecAlgo(bool correctForTimeslew, bool correctForPulse, float phaseNS) : 
  correctForTimeslew_(correctForTimeslew),
  correctForPulse_(correctForPulse),
  phaseNS_(phaseNS), runnum_(0), setLeakCorrection_(false), puCorrMethod_(0)
{ 
  
  pulseCorr_ = std::auto_ptr<HcalPulseContainmentManager>(
	       new HcalPulseContainmentManager(MaximumFractionalError)
							  );
}
  

HcalSimpleRecAlgo::HcalSimpleRecAlgo() : 
  correctForTimeslew_(false), runnum_(0), puCorrMethod_(0) { }


void HcalSimpleRecAlgo::beginRun(edm::EventSetup const & es)
{
  pulseCorr_->beginRun(es);
}


void HcalSimpleRecAlgo::endRun()
{
  pulseCorr_->endRun();
}


void HcalSimpleRecAlgo::initPulseCorr(int toadd) {
}

void HcalSimpleRecAlgo::setRecoParams(bool correctForTimeslew, bool correctForPulse, bool setLeakCorrection, int pileupCleaningID, float phaseNS){
   correctForTimeslew_=correctForTimeslew;
   correctForPulse_=correctForPulse;
   phaseNS_=phaseNS;
   setLeakCorrection_=setLeakCorrection;
   pileupCleaningID_=pileupCleaningID;
}

void HcalSimpleRecAlgo::setForData (int runnum) { runnum_ = runnum;}

void HcalSimpleRecAlgo::setLeakCorrection () { setLeakCorrection_ = true;}

void HcalSimpleRecAlgo::setHBHEPileupCorrection(
     boost::shared_ptr<AbsOOTPileupCorrection> corr)
{
    hbhePileupCorr_ = corr;
}

void HcalSimpleRecAlgo::setHFPileupCorrection(
     boost::shared_ptr<AbsOOTPileupCorrection> corr)
{
    hfPileupCorr_ = corr;
}

void HcalSimpleRecAlgo::setHOPileupCorrection(
     boost::shared_ptr<AbsOOTPileupCorrection> corr)
{
    hoPileupCorr_ = corr;
}

void HcalSimpleRecAlgo::setBXInfo(const BunchXParameter* info,
                                  const unsigned lenInfo)
{
    bunchCrossingInfo_ = info;
    lenBunchCrossingInfo_ = lenInfo;
}

///Timeshift correction for HPDs based on the position of the peak ADC measurement.
///  Allows for an accurate determination of the relative phase of the pulse shape from
///  the HPD.  Calculated based on a weighted sum of the -1,0,+1 samples relative to the peak
///  as follows:  wpksamp = (0*sample[0] + 1*sample[1] + 2*sample[2]) / (sample[0] + sample[1] + sample[2])
///  where sample[1] is the maximum ADC sample value.
static float timeshift_ns_hbheho(float wpksamp);

///Same as above, but for the HF PMTs.
static float timeshift_ns_hf(float wpksamp);

/// Ugly hack to apply energy corrections to some HB- cells
static float eCorr(int ieta, int iphi, double ampl, int runnum);

/// Leak correction 
static float leakCorr(double energy);


namespace HcalSimpleRecAlgoImpl {
  template<class Digi>
  inline float recoHFTime(const Digi& digi, const int maxI, const double amp_fC,
                          const bool slewCorrect, double maxA, float t0, float t2)
  {
    // Handle negative excursions by moving "zero":
    float zerocorr=std::min(t0,t2);
    if (zerocorr<0.f) {
      t0   -= zerocorr;
      t2   -= zerocorr;
      maxA -= zerocorr;
    }
    
    // pair the peak with the larger of the two neighboring time samples
    float wpksamp=0.f;
    if (t0>t2) {
      wpksamp = t0+maxA;
      if (wpksamp != 0.f) wpksamp = maxA/wpksamp;
    } else {
      wpksamp = maxA+t2;
      if (wpksamp != 0.f) wpksamp = 1.+(t2/wpksamp);
    }

    float time = (maxI - digi.presamples())*25.0 + timeshift_ns_hf(wpksamp);

    if (slewCorrect && amp_fC > 0.0) {
      // -5.12327 - put in calibs.timecorr()
      double tslew=exp(0.337681-5.94689e-4*amp_fC)+exp(2.44628-1.34888e-2*amp_fC);
      time -= (float)tslew;
    }

    return time;
  }


  template<class Digi>
  inline void removePileup(const Digi& digi, const HcalCoder& coder,
                           const HcalCalibrations& calibs,
                           const int ifirst, const int n,
                           const bool pulseCorrect,
                           const HcalPulseContainmentCorrection* corr,
                           const AbsOOTPileupCorrection* pileupCorrection,
                           const BunchXParameter* bxInfo, const unsigned lenInfo,
                           double* p_maxA, double* p_ampl, double* p_uncorr_ampl,
                           double* p_fc_ampl, int* p_nRead, int* p_maxI,
                           bool* leakCorrApplied, float* p_t0, float* p_t2)
  {
    CaloSamples cs;
    coder.adc2fC(digi,cs);
    const int nRead = cs.size();
    const int iStop = std::min(nRead, n + ifirst);

    // Signal energy will be calculated both with
    // and without OOT pileup corrections. Try to
    // arrange the calculations so that we do not
    // repeat them.
    double uncorrectedEnergy[CaloSamples::MAXSAMPLES], buf[CaloSamples::MAXSAMPLES];
    double* correctedEnergy = 0;
    double fc_ampl = 0.0, corr_fc_ampl = 0.0;
    bool pulseShapeCorrApplied = false, readjustTiming = false;
    *leakCorrApplied = false;

    if (pileupCorrection)
    {
        correctedEnergy = &buf[0];

        double correctionInput[CaloSamples::MAXSAMPLES];
        double gains[CaloSamples::MAXSAMPLES];

        for (int i=0; i<nRead; ++i)
        {
            const int capid = digi[i].capid();
            correctionInput[i] = cs[i] - calibs.pedestal(capid);
            gains[i] = calibs.respcorrgain(capid);
        }

        for (int i=ifirst; i<iStop; ++i)
            fc_ampl += correctionInput[i];

        const bool useGain = pileupCorrection->inputIsEnergy();
        for (int i=0; i<nRead; ++i)
        {
            uncorrectedEnergy[i] = correctionInput[i]*gains[i];
            if (useGain)
                correctionInput[i] = uncorrectedEnergy[i];
        }

        pileupCorrection->apply(digi.id(), correctionInput, nRead,
                                bxInfo, lenInfo, ifirst, n,
                                correctedEnergy, CaloSamples::MAXSAMPLES,
                                &pulseShapeCorrApplied, leakCorrApplied,
                                &readjustTiming);
        if (useGain)
        {
            // Gain factors have been already applied.
            // Divide by them for accumulating corr_fc_ampl.
            for (int i=ifirst; i<iStop; ++i)
                if (gains[i])
                    corr_fc_ampl += correctedEnergy[i]/gains[i];
        }
        else
        {
            for (int i=ifirst; i<iStop; ++i)
                corr_fc_ampl += correctedEnergy[i];
            for (int i=0; i<nRead; ++i)
                correctedEnergy[i] *= gains[i];
        }
    }
    else
    {
        correctedEnergy = &uncorrectedEnergy[0];

        // In this situation, we do not need to process all time slices
        const int istart = std::max(ifirst - 1, 0);
        const int iend = std::min(n + ifirst + 1, nRead);
        for (int i=istart; i<iend; ++i)
        {
            const int capid = digi[i].capid();
            float ta = cs[i] - calibs.pedestal(capid);
            if (i >= ifirst && i < iStop)
                fc_ampl += ta;
            ta *= calibs.respcorrgain(capid);
            uncorrectedEnergy[i] = ta;
        }
        corr_fc_ampl = fc_ampl;
    }

    // Uncorrected and corrected energies
    double ampl = 0.0, corr_ampl = 0.0;
    for (int i=ifirst; i<iStop; ++i)
    {
        ampl += uncorrectedEnergy[i];
        corr_ampl += correctedEnergy[i];
    }

    // Apply phase-based amplitude correction:
    if (corr && pulseCorrect)
    {
        ampl *= corr->getCorrection(fc_ampl);
        if (pileupCorrection)
        {
            if (!pulseShapeCorrApplied)
                corr_ampl *= corr->getCorrection(corr_fc_ampl);
        }
        else
            corr_ampl = ampl;
    }

    // Which energies we want to use for timing?
    const double *etime = readjustTiming ? &correctedEnergy[0] : &uncorrectedEnergy[0];
    int maxI = -1; double maxA = -1.e300;
    for (int i=ifirst; i<iStop; ++i)
        if (etime[i] > maxA)
        {
            maxA = etime[i];
            maxI = i;
        }

    // Fill out the output
    *p_maxA = maxA;
    *p_ampl = corr_ampl;
    *p_uncorr_ampl = ampl;
    *p_fc_ampl = readjustTiming ? corr_fc_ampl : fc_ampl;
    *p_nRead = nRead;
    *p_maxI = maxI;

    if (maxI <= 0 || maxI >= (nRead-1))
    {
      LogDebug("HCAL Pulse") << "HcalSimpleRecAlgoImpl::removePileup :" 
					       << " Invalid max amplitude position, " 
					       << " max Amplitude: " << maxI
					       << " first: " << ifirst
					       << " last: " << ifirst + n
					       << std::endl;
      *p_t0 = 0.f;
      *p_t2 = 0.f;
    }
    else
    {
      *p_t0 = etime[maxI - 1];
      *p_t2 = etime[maxI + 1];
    }
  }

//  HcalPulseShapes::Shape hpdshape;
//  double fitFunction(double* x, double* pars);
  double funcHPDShapeData(double* x, double* pars);
  double funcHPDShapeMC(double* x, double* pars);
  double func_DoublePulse_HPDShapeData(double* x, double* pars);
  double func_DoublePulse_HPDShapeMC(double* x, double* pars);
  void fcn1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pars, Int_t iflag);
  void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pars, Int_t iflag);
// Note that: charge is without pedestal subtraction!
  int pulseShapeFit(const std::vector<double> &charge, const std::vector<double> &ped, const double TSTOT, std::vector<double> &fitParsVec);
  double psFit_x[10], psFit_y[10], psFit_erry[10];
  int psFit_isData;

  template<class Digi, class RecHit>
  inline RecHit reco(const Digi& digi, const HcalCoder& coder,
                     const HcalCalibrations& calibs, 
		     const int ifirst, const int n, const bool slewCorrect,
                     const bool pulseCorrect, const HcalPulseContainmentCorrection* corr,
		     const HcalTimeSlew::BiasSetting slewFlavor,
                     const int runnum, const bool useLeak,
                     const AbsOOTPileupCorrection* pileupCorrection,
                     const BunchXParameter* bxInfo, const unsigned lenInfo, const int puCorrMethod)
  {
    double fc_ampl, ampl, uncorr_ampl, maxA;
    int nRead, maxI;
    bool leakCorrApplied;
    float t0, t2;

    removePileup(digi, coder, calibs, ifirst, n,
                 pulseCorrect, corr, pileupCorrection,
                 bxInfo, lenInfo, &maxA, &ampl,
                 &uncorr_ampl, &fc_ampl, &nRead, &maxI,
                 &leakCorrApplied, &t0, &t2);

    float time = -9999;
    if (maxI > 0 && maxI < (nRead - 1))
    {
      // Handle negative excursions by moving "zero":
      float minA=t0;
      if (maxA<minA) minA=maxA;
      if (t2<minA)   minA=t2;
      if (minA<0) { maxA-=minA; t0-=minA; t2-=minA; } // positivizes all samples

      float wpksamp = (t0 + maxA + t2);
      if (wpksamp!=0) wpksamp=(maxA + 2.0*t2) / wpksamp; 
      time = (maxI - digi.presamples())*25.0 + timeshift_ns_hbheho(wpksamp);

      if (slewCorrect) time-=HcalTimeSlew::delay(std::max(1.0,fc_ampl),slewFlavor);

      time=time-calibs.timecorr(); // time calibration
    }

    // Temporary hack to apply energy-dependent corrections to some HB- cells
    if (runnum > 0) {
      const HcalDetId& cell = digi.id();
      if (cell.subdet() == HcalBarrel) {
        const int ieta = cell.ieta();
        const int iphi = cell.iphi();
        ampl *= eCorr(ieta, iphi, ampl, runnum);
        uncorr_ampl *= eCorr(ieta, iphi, uncorr_ampl, runnum);
      }
    }

    // Correction for a leak to pre-sample
    if(useLeak && !leakCorrApplied) {
      ampl *= leakCorr(ampl); 
      uncorr_ampl *= leakCorr(uncorr_ampl); 
    }

    if( puCorrMethod == 2 ){
       if( runnum >0 /*data*/ ){
          psFit_isData = 1;
       }else{
          psFit_isData = 0;
       }

       CaloSamples cs;
       coder.adc2fC(digi,cs);
       std::vector<double> charge, ped;
       double TSTOT = 0;
       for(int ip=0; ip<cs.size(); ip++){
          charge.push_back(cs[ip]);

          const int capid = digi[ip].capid();
          double perped = calibs.pedestal(capid);
          ped.push_back(perped);

          TSTOT += cs[ip] - perped;
       }
       std::vector<double> fitParsVec;
       pulseShapeFit(charge, ped, TSTOT, fitParsVec);
       time = fitParsVec[1]; ampl = fitParsVec[0]; uncorr_ampl = fitParsVec[0];
    }

    RecHit rh(digi.id(),ampl,time);
    setRawEnergy(rh, static_cast<float>(uncorr_ampl));
    return rh;
  }

  int pulseShapeFit(const std::vector<double> &charge, const std::vector<double> &ped, const double TSTOT, std::vector<double> &fitParsVec){
  
     int n_max=0;
     int n_above_thr=0;
     int first_above_thr_index=-1;
     int max_index[10]={0,0,0,0,0,0,0,0,0,0};
 
     double TSMAX=0;
     double TSMAX_NOPED=0;
     int i_tsmax=0;
 
     for(int i=0;i<10;i++){
        if(charge[i]>TSMAX){
           TSMAX=charge[i];
           TSMAX_NOPED=charge[i]-ped[i];
           i_tsmax = i;
        }
     }
 
     double TIMES[10]={-100,-75,-50,-25,0,25,50,75,100,125};
 
     if(n_max==0){
        max_index[0]=i_tsmax;
     }
    
     double error = 1.;
     for(int i=0;i<10;i++){
        psFit_x[i]=i;
        psFit_y[i]=charge[i];
        psFit_erry[i]=error;
     }
    
     TMinuit * gMinuit = new TMinuit(5);
     gMinuit->SetPrintLevel(-1);
  
     for(int i=0;i!=10;++i){
        if((charge[i])>6){
           n_above_thr++;
           if(first_above_thr_index==-1){
              first_above_thr_index=i;
           }
        }
     }
     for(int i=1;i!=9;++i){
        if(charge[i-1]>6 && charge[i]>6 && charge[i+1]>6){
           if(charge[i-1]<charge[i] && charge[i+1]<charge[i]){
              max_index[n_max]=i;
              n_max++;
           }
        }
     }
     if(n_max==0){
        max_index[0]=i_tsmax;
     }
 
     if(n_above_thr<=5){
        gMinuit->SetFCN(fcn1);
        double arglist[10];
        int ierflg = 0;
 
        arglist[0] = 1;
        gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
 
        // Set starting values and step sizes for parameters
        double vstart[3] = {TIMES[i_tsmax-1],TSMAX_NOPED,0};
        double step[3] = {0.1,0.1,0.1};
        gMinuit->mnparm(0, "time", vstart[0], step[0], -100,75,ierflg);
        gMinuit->mnparm(1, "energy", vstart[1], step[1], 0,TSTOT,ierflg);
        gMinuit->mnparm(2, "ped", vstart[2], step[2], 0,0,ierflg);
 
        double chi2=9999.;
 
        for(int tries=0; tries<=3;tries++){
           // Now ready for minimization step
           arglist[0] = 500;
           arglist[1] = 1.;
           gMinuit->mnexcm("Migrad", arglist ,2,ierflg);
 
           double chi2valfit,edm,errdef;
           int nvpar,nparx,icstat;
           gMinuit->mnstat(chi2valfit,edm,errdef,nvpar,nparx,icstat);
 
           if(chi2>chi2valfit+0.01){
              chi2=chi2valfit;
              if(tries==0){
                 arglist[0]=0;
                 gMinuit->mnexcm("SCAN",arglist,1,ierflg);
              } else if(tries==1){
                 arglist[0] = 1;
                 gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
              } else if(tries==2){
                 arglist[0] = 2;
                 gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
              }
           } else {
              break;
           }
        }
     }else{
        gMinuit->SetFCN(fcn2);
        double arglist[10];
        int ierflg = 0;
 
        arglist[0] = 1;
        gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
 
        if(n_max==1){
        // Set starting values and step sizes for parameters
           double vstart[5] = {TIMES[i_tsmax-1],TIMES[first_above_thr_index-1],TSMAX_NOPED,0,0};
           Double_t step[5] = {0.1,0.1,0.1,0.1,0.1};
           gMinuit->mnparm(0, "time1", vstart[0], step[0], -100,75,ierflg);
           gMinuit->mnparm(1, "time2", vstart[1], step[1], -100,75,ierflg);
           gMinuit->mnparm(2, "energy1", vstart[2], step[2], 0,TSTOT,ierflg);
           gMinuit->mnparm(3, "energy2", vstart[3], step[3], 0,TSTOT,ierflg);
           gMinuit->mnparm(4, "ped", vstart[4], step[4], 0,0,ierflg);
 
           double chi2=9999.;
 
           for(int tries=0; tries<=3;tries++) {
           // Now ready for minimization step
              arglist[0] = 500;
              arglist[1] = 1.;
              gMinuit->mnexcm("Migrad", arglist ,2,ierflg);
 
              double chi2valfit,edm,errdef;
              int nvpar,nparx,icstat;
              gMinuit->mnstat(chi2valfit,edm,errdef,nvpar,nparx,icstat);
 
              if(chi2>chi2valfit+0.01) {
                 chi2=chi2valfit;
                 if(tries==0) {
                    arglist[0]=0;
                    gMinuit->mnexcm("SCAN",arglist,1,ierflg);
                 } else if(tries==1) {
                    arglist[0] = 1;
                    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
                 } else if(tries==2) {
                    arglist[0] = 2;
                    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
                 }
              } else {
                 break;
              }
           }
        } else if(n_max>=2) {
        // Set starting values and step sizes for parameters
           double vstart[5] = {TIMES[max_index[0]-1],TIMES[max_index[1]-1],TSMAX_NOPED,0,0};
           double step[5] = {0.1,0.1,0.1,0.1,0.1};
           gMinuit->mnparm(0, "time1", vstart[0], step[0], -100,75,ierflg);
           gMinuit->mnparm(1, "time2", vstart[1], step[1], -100,75,ierflg);
           gMinuit->mnparm(2, "energy1", vstart[2], step[2], 0,TSTOT,ierflg);
           gMinuit->mnparm(3, "energy2", vstart[3], step[3], 0,TSTOT,ierflg);
           gMinuit->mnparm(4, "ped", vstart[4], step[4], 0,0,ierflg);
 
           double chi2=9999.;
 
           for(int tries=0; tries<=3;tries++) {
           // Now ready for minimization step
              arglist[0] = 500;
              arglist[1] = 1.;
              gMinuit->mnexcm("Migrad", arglist ,2,ierflg);
 
              double chi2valfit,edm,errdef;
              int nvpar,nparx,icstat;
              gMinuit->mnstat(chi2valfit,edm,errdef,nvpar,nparx,icstat);
 
              if(chi2>chi2valfit+0.01) {
                 chi2=chi2valfit;
                 if(tries==0) {
                    arglist[0]=0;
                    gMinuit->mnexcm("SCAN",arglist,1,ierflg);
                 } else if(tries==1) {
                    arglist[0] = 1;
                    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
                 } else if(tries==2) {
                    arglist[0] = 2;
                    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
                 }
              } else {
                 break;
              }
           }
        }
     }
 
     int fitStatus=gMinuit->GetStatus();
 
     double timeval1fit=-999;
     double timeval1fit_err=-999;
     double chargeval1fit=-999;
     double chargeval1fit_err=-999;
     double timeval2fit=-999;
     double timeval2fit_err=-999;
     double chargeval2fit=-999;
     double chargeval2fit_err=-999;
     double pedvalfit=-999;
     double pedvalfit_err=-999;

     if(n_above_thr<=5) {
        gMinuit->GetParameter(0,timeval1fit,timeval1fit_err);
        gMinuit->GetParameter(1,chargeval1fit,chargeval1fit_err);
        gMinuit->GetParameter(2,pedvalfit,pedvalfit_err);
     } else {
        gMinuit->GetParameter(0,timeval1fit,timeval1fit_err);
        gMinuit->GetParameter(1,timeval2fit,timeval2fit_err);
        gMinuit->GetParameter(2,chargeval1fit,chargeval1fit_err);
        gMinuit->GetParameter(3,chargeval2fit,chargeval2fit_err);
        gMinuit->GetParameter(4,pedvalfit,pedvalfit_err);
     }

     double chi2valfit,edm,errdef;
     int nvpar,nparx,icstat;
     gMinuit->mnstat(chi2valfit,edm,errdef,nvpar,nparx,icstat);

     double timevalfit=0.;
     double chargevalfit=0.;
     if(n_above_thr<=5) {
        timevalfit=timeval1fit;
        chargevalfit=chargeval1fit;
     } else { 
       // if timeval1fit and timeval2fit are differnt, choose the one which is closer to zero
        if(abs(timeval1fit)<abs(timeval2fit)) {
           timevalfit=timeval1fit;
           chargevalfit=chargeval1fit;
        } else if(abs(timeval2fit)<abs(timeval1fit)) {
           timevalfit=timeval2fit;
           chargevalfit=chargeval2fit;
        }
        // if the two times are the same, then for charge we just sum the two
        else if(timeval1fit==timeval2fit) {
           timevalfit=(timeval1fit+timeval2fit)/2;
           chargevalfit=chargeval1fit+chargeval2fit;
        } else {
           timevalfit=-999.;
           chargevalfit=-999.;
        }
     }

     if( gMinuit ) delete gMinuit;

     fitParsVec.clear();

     fitParsVec.push_back(chargevalfit);
     fitParsVec.push_back(timevalfit);
     fitParsVec.push_back(pedvalfit);
     fitParsVec.push_back(chi2valfit);

     return fitStatus;
  }
 
  void fcn1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pars, Int_t iflag) {
     const int nbins = 10;
     int i;

     //calculate chisquare
     double chisq = 0;
     double delta;
     double val[1];
     for (i=0;i<nbins; i++) {
        val[0]=psFit_x[i];
        delta = (psFit_y[i]- funcHPDShapeData(val,pars))/psFit_erry[i];
        if( !psFit_isData ){
           delta = (psFit_y[i]- funcHPDShapeMC(val,pars))/psFit_erry[i];
        }
        chisq += delta*delta;
     }
     f = chisq;
  }

  void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pars, Int_t iflag) {
     const int nbins = 10;
     int i;

     //calculate chisquare
     double chisq = 0;
     double delta;
     double val[1];
     for (i=0;i<nbins; i++) {
        val[0] = psFit_x[i];
        delta = (psFit_y[i] - func_DoublePulse_HPDShapeData(val,pars))/psFit_erry[i];
        if( !psFit_isData ){
           delta = (psFit_y[i] - func_DoublePulse_HPDShapeMC(val,pars))/psFit_erry[i];
        }
        chisq += delta*delta;
     }
     f = chisq;
  }

  double funcHPDShapeData(double* x, double* pars) {

  // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
       int nbin = 256;
     
     std::vector<float> shape(nbin,0.0);
     
     shape[0]=0;
     shape[1]=0.00033394;
     shape[2]=0.00136382;
     shape[3]=0.00345003;
     shape[4]=0.00692451;
     shape[5]=0.0114013;
     shape[6]=0.0164575;
     shape[7]=0.0216718;
     shape[8]=0.0266611;
     shape[9]=0.0311094;
     shape[10]=0.0347871;
     shape[11]=0.0375579;
     shape[12]=0.0393753;
     shape[13]=0.0402696;
     shape[14]=0.0403301;
     shape[15]=0.0396845;
     shape[16]=0.0384793;
     shape[17]=0.0368628;
     shape[18]=0.0349727;
     shape[19]=0.032928;
     shape[20]=0.0308253;
     shape[21]=0.0287377;
     shape[22]=0.0267174;
     shape[23]=0.024798;
     shape[24]=0.022999;
     shape[25]=0.0213288;
     shape[26]=0.0197884;
     shape[27]=0.0183738;
     shape[28]=0.0170779;
     shape[29]=0.0158922;
     shape[30]=0.0148077;
     shape[31]=0.0138152;
     shape[32]=0.0129064;
     shape[33]=0.012073;
     shape[34]=0.011308;
     shape[35]=0.0106047;
     shape[36]=0.00995714;
     shape[37]=0.00936001;
     shape[38]=0.00880855;
     shape[39]=0.00829848;
     shape[40]=0.007826;
     shape[41]=0.00738768;
     shape[42]=0.00698044;
     shape[43]=0.00660155;
     shape[44]=0.00624853;
     shape[45]=0.00591916;
     shape[46]=0.00561146;
     shape[47]=0.00532362;
     shape[48]=0.00505403;
     shape[49]=0.00480122;
     shape[50]=0.00456388;
     shape[51]=0.00434079;
     shape[52]=0.0041309;
     shape[53]=0.0039332;
     shape[54]=0.00374682;
     shape[55]=0.00357093;
     shape[56]=0.00340479;
     shape[57]=0.00324773;
     shape[58]=0.00309914;
     shape[59]=0.00295845;
     shape[60]=0.00282513;
     shape[61]=0.00269872;
     shape[62]=0.00257878;
     shape[63]=0.0024649;
     shape[64]=0.00235672;
     shape[65]=0.00225389;
     shape[66]=0.0021561;
     shape[67]=0.00206304;
     shape[68]=0.00197446;
     shape[69]=0.00189009;
     shape[70]=0.00180971;
     shape[71]=0.00173309;
     shape[72]=0.00166003;
     shape[73]=0.00159034;
     shape[74]=0.00152384;
     shape[75]=0.00146036;
     shape[76]=0.00139975;
     shape[77]=0.00134186;
     shape[78]=0.00128656;
     shape[79]=0.00123371;
     shape[80]=0.0011832;
     shape[81]=0.0011349;
     shape[82]=0.00108872;
     shape[83]=0.00104454;
     shape[84]=0.00100228;
     shape[85]=0.000961835;
     shape[86]=0.00092313;
     shape[87]=0.000886081;
     shape[88]=0.000850609;
     shape[89]=0.000816643;
     shape[90]=0.000784113;
     shape[91]=0.000752953;
     shape[92]=0.000723102;
     shape[93]=0.0006945;
     shape[94]=0.000667091;
     shape[95]=0.000640823;
     shape[96]=0.000615643;
     shape[97]=0.000591505;
     shape[98]=0.000568362;
     shape[99]=0.00054617;
     shape[100]=0.000524888;
     shape[101]=0.000504476;
     shape[102]=0.000484897;
     shape[103]=0.000466114;
     shape[104]=0.000448094;
     shape[105]=0.000430803;
     shape[106]=0.000414211;
     shape[107]=0.000398286;
     shape[108]=0.000383002;
     shape[109]=0.000368331;
     shape[110]=0.000354248;
     shape[111]=0.000340726;
     shape[112]=0.000327743;
     shape[113]=0.000315276;
     shape[114]=0.000303304;
     shape[115]=0.000291806;
     shape[116]=0.000280762;
     shape[117]=0.000270153;
     shape[118]=0.000259962;
     shape[119]=0.000250171;
     shape[120]=0.000240764;
     shape[121]=0.000231725;
     shape[122]=0.000223038;
     shape[123]=0.000214691;
     shape[124]=0.000206667;
     shape[125]=0.000198956;
     shape[126]=0.000191543;
     shape[127]=0.000184416;
     shape[128]=0.000177565;
     shape[129]=0.000170978;
     shape[130]=0.000164644;
     shape[131]=0.000158554;
     shape[132]=0.000152697;
     shape[133]=0.000147064;
     shape[134]=0.000141646;
     shape[135]=0.000136435;
     shape[136]=0.000131423;
     shape[137]=0.000126601;
     shape[138]=0.000121962;
     shape[139]=0.000117498;
     shape[140]=0.000113204;
     shape[141]=0.000109071;
     shape[142]=0.000105095;
     shape[143]=0.000101268;
     shape[144]=9.75848e-05;
     shape[145]=9.404e-05;
     shape[146]=9.0628e-05;
     shape[147]=8.73436e-05;
     shape[148]=8.4182e-05;
     shape[149]=8.11383e-05;
     shape[150]=7.8208e-05;
     shape[151]=7.53866e-05;
     shape[152]=7.26701e-05;
     shape[153]=7.00543e-05;
     shape[154]=6.75353e-05;
     shape[155]=6.51096e-05;
     shape[156]=6.27734e-05;
     shape[157]=6.05233e-05;
     shape[158]=5.83562e-05;
     shape[159]=5.62687e-05;
     shape[160]=5.42579e-05;
     shape[161]=5.23209e-05;
     shape[162]=5.04549e-05;
     shape[163]=4.86571e-05;
     shape[164]=4.69251e-05;
     shape[165]=4.52562e-05;
     shape[166]=4.36482e-05;
     shape[167]=4.20987e-05;
     shape[168]=4.06056e-05;
     shape[169]=3.91667e-05;
     shape[170]=3.778e-05;
     shape[171]=3.64436e-05;
     shape[172]=3.51555e-05;
     shape[173]=3.3914e-05;
     shape[174]=3.27173e-05;
     shape[175]=3.14624e-05;
     shape[176]=3.00317e-05;
     shape[177]=2.83024e-05;
     shape[178]=2.61543e-05;
     shape[179]=2.36806e-05;
     shape[180]=2.09939e-05;
     shape[181]=1.82145e-05;
     shape[182]=1.5459e-05;
     shape[183]=1.28302e-05;
     shape[184]=1.041e-05;
     shape[185]=8.25538e-06;
     shape[186]=6.3974e-06;
     shape[187]=4.84374e-06;
     shape[188]=3.58271e-06;
     shape[189]=2.58848e-06;
     shape[190]=1.82658e-06;
     shape[191]=1.25879e-06;
     shape[192]=8.47151e-07;
     shape[193]=5.56713e-07;
     shape[194]=3.57223e-07;
     shape[195]=2.23802e-07;
     shape[196]=1.36893e-07;
     shape[197]=8.17475e-08;
     shape[198]=4.76566e-08;
     shape[199]=2.71208e-08;
     shape[200]=1.50656e-08;
     shape[201]=8.16834e-09;
     shape[202]=4.32189e-09;
     shape[203]=2.2309e-09;
     shape[204]=1.12281e-09;
     shape[205]=5.50335e-10;
     shape[206]=2.62018e-10;
     shape[207]=1.20477e-10;
     shape[208]=5.27642e-11;
     shape[209]=2.12197e-11;
     shape[210]=7.34902e-12;
     shape[211]=1.76216e-12;
     shape[212]=0;
     shape[213]=0;
     shape[214]=0;
     shape[215]=0;
     shape[216]=0;
     shape[217]=0;
     shape[218]=0;
     shape[219]=0;
     shape[220]=0;
     shape[221]=0;
     shape[222]=0;
     shape[223]=0;
     shape[224]=0;
     shape[225]=0;
     shape[226]=0;
     shape[227]=0;
     shape[228]=0;
     shape[229]=0;
     shape[230]=0;
     shape[231]=0;
     shape[232]=0;
     shape[233]=0;
     shape[234]=0;
     shape[235]=0;
     shape[236]=0;
     shape[237]=0;
     shape[238]=0;
     shape[239]=0;
     shape[240]=0;
     shape[241]=0;
     shape[242]=0;
     shape[243]=0;
     shape[244]=0;
     shape[245]=0;
     shape[246]=0;
     shape[247]=0;
     shape[248]=0;
     shape[249]=0;
     shape[250]=0;
     shape[251]=0;
     shape[252]=0;
     shape[253]=0;
     shape[254]=0;
     shape[255]=0;
     
     std::vector<float> ntmpbin(10,0.0);  // zeroing output binned pulse shape
     std::vector<float> ntmpshift(nbin,0.0);  // zeroing output shifted pulse shape  
     
     int xx = (int)x[0];
     double w1 = pars[0];
     double w2 = pars[1];
     double w3 = pars[2];
        
     TH1F *h1=new TH1F("h1","test",256,0,256);
     
     for(int i=0;i<nbin;i++){
        h1->SetBinContent(i+1,shape[i]);
     }
     
     for(int i=0;i<nbin;i++){
        if(i<w1+98.5) {
     	   ntmpshift[i]=0;
        } else {
     	   ntmpshift[i]=h1->Interpolate(i-98.5-w1);
        }
    }
     
    for(int i=0;i<nbin;i++) {
       if(i<250) {
     	  if(i<25) {
     	     ntmpbin[0]+=ntmpshift[i];
     	  } else if(i>=25&&i<50) {
     	     ntmpbin[1]+=ntmpshift[i];
     	  } else if(i>=50&&i<75) {
     	     ntmpbin[2]+=ntmpshift[i];
     	  } else if(i>=75&&i<100) {
     	     ntmpbin[3]+=ntmpshift[i];
     	  } else if(i>=100&&i<125) {
     	     ntmpbin[4]+=ntmpshift[i];
     	  } else if(i>=125&&i<150) {
     	     ntmpbin[5]+=ntmpshift[i];
     	  } else if(i>=150&&i<175) {
     	     ntmpbin[6]+=ntmpshift[i];
     	  } else if(i>=175&&i<200) {
     	     ntmpbin[7]+=ntmpshift[i];
     	  } else if(i>=200&&i<225) {
     	     ntmpbin[8]+=ntmpshift[i];
     	  } else if(i>=225&&i<250) {
     	     ntmpbin[9]+=ntmpshift[i];
     	  }
       }
    }
    
    if( h1 ) delete h1; 
      
    return w2*ntmpbin[xx]+w3;
  }

  double funcHPDShapeMC(double* x, double* pars) {
     
       // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
     int nbin = 256;
     
     
     std::vector<float> shape(nbin,0.0);
     
     shape[0]=0;
     shape[1]=0.000572922;
     shape[2]=0.00231338;
     shape[3]=0.00576442;
     shape[4]=0.0113605;
     shape[5]=0.0182484;
     shape[6]=0.025526;
     shape[7]=0.0323732;
     shape[8]=0.0381537;
     shape[9]=0.042472;
     shape[10]=0.0451786;
     shape[11]=0.0463339;
     shape[12]=0.0461474;
     shape[13]=0.0449083;
     shape[14]=0.0429253;
     shape[15]=0.0404816;
     shape[16]=0.0378101;
     shape[17]=0.0350845;
     shape[18]=0.0324225;
     shape[19]=0.0298955;
     shape[20]=0.0275404;
     shape[21]=0.0253707;
     shape[22]=0.0233854;
     shape[23]=0.0215755;
     shape[24]=0.0199283;
     shape[25]=0.0184295;
     shape[26]=0.0170654;
     shape[27]=0.0158227;
     shape[28]=0.0146895;
     shape[29]=0.0136548;
     shape[30]=0.012709;
     shape[31]=0.0118431;
     shape[32]=0.0110495;
     shape[33]=0.0103211;
     shape[34]=0.00965171;
     shape[35]=0.00903563;
     shape[36]=0.00846789;
     shape[37]=0.00794399;
     shape[38]=0.00745989;
     shape[39]=0.00701198;
     shape[40]=0.00659701;
     shape[41]=0.00621204;
     shape[42]=0.00585447;
     shape[43]=0.00552192;
     shape[44]=0.00521226;
     shape[45]=0.00492357;
     shape[46]=0.00465411;
     shape[47]=0.00440231;
     shape[48]=0.00416676;
     shape[49]=0.00394616;
     shape[50]=0.00373936;
     shape[51]=0.0035453;
     shape[52]=0.00336301;
     shape[53]=0.00319162;
     shape[54]=0.00303033;
     shape[55]=0.00287843;
     shape[56]=0.00273523;
     shape[57]=0.00260015;
     shape[58]=0.00247262;
     shape[59]=0.00235214;
     shape[60]=0.00223823;
     shape[61]=0.00213048;
     shape[62]=0.00202847;
     shape[63]=0.00193186;
     shape[64]=0.00184031;
     shape[65]=0.0017535;
     shape[66]=0.00167114;
     shape[67]=0.00159298;
     shape[68]=0.00151877;
     shape[69]=0.00144827;
     shape[70]=0.00138128;
     shape[71]=0.00131759;
     shape[72]=0.00125703;
     shape[73]=0.00119941;
     shape[74]=0.00114459;
     shape[75]=0.0010924;
     shape[76]=0.00104271;
     shape[77]=0.000995392;
     shape[78]=0.000950312;
     shape[79]=0.00090736;
     shape[80]=0.000866425;
     shape[81]=0.000827405;
     shape[82]=0.000790203;
     shape[83]=0.000754728;
     shape[84]=0.000720894;
     shape[85]=0.000688621;
     shape[86]=0.000657831;
     shape[87]=0.000628453;
     shape[88]=0.000600418;
     shape[89]=0.000573662;
     shape[90]=0.000548124;
     shape[91]=0.000523744;
     shape[92]=0.000500469;
     shape[93]=0.000478247;
     shape[94]=0.000457027;
     shape[95]=0.000436763;
     shape[96]=0.000417411;
     shape[97]=0.000398927;
     shape[98]=0.000381273;
     shape[99]=0.000364409;
     shape[100]=0.000348299;
     shape[101]=0.000332909;
     shape[102]=0.000318206;
     shape[103]=0.000304158;
     shape[104]=0.000290736;
     shape[105]=0.00027791;
     shape[106]=0.000265655;
     shape[107]=0.000253945;
     shape[108]=0.000242753;
     shape[109]=0.000232059;
     shape[110]=0.000221838;
     shape[111]=0.00021207;
     shape[112]=0.000202734;
     shape[113]=0.000193811;
     shape[114]=0.000185283;
     shape[115]=0.000177132;
     shape[116]=0.000169341;
     shape[117]=0.000161894;
     shape[118]=0.000154775;
     shape[119]=0.000147971;
     shape[120]=0.000141466;
     shape[121]=0.000135249;
     shape[122]=0.000129305;
     shape[123]=0.000123624;
     shape[124]=0.000118192;
     shape[125]=0.000113;
     shape[126]=0.000108037;
     shape[127]=0.000103292;
     shape[128]=9.87554e-05;
     shape[129]=9.44188e-05;
     shape[130]=9.02728e-05;
     shape[131]=8.63093e-05;
     shape[132]=8.252e-05;
     shape[133]=7.8481e-05;
     shape[134]=7.37323e-05;
     shape[135]=6.78316e-05;
     shape[136]=6.03989e-05;
     shape[137]=5.19781e-05;
     shape[138]=4.31841e-05;
     shape[139]=3.46072e-05;
     shape[140]=2.67337e-05;
     shape[141]=1.98962e-05;
     shape[142]=1.42599e-05;
     shape[143]=9.83903e-06;
     shape[144]=6.53358e-06;
     shape[145]=4.17456e-06;
     shape[146]=2.56593e-06;
     shape[147]=1.51697e-06;
     shape[148]=8.62467e-07;
     shape[149]=4.715e-07;
     shape[150]=2.47822e-07;
     shape[151]=1.25217e-07;
     shape[152]=6.08124e-08;
     shape[153]=2.83828e-08;
     shape[154]=1.27271e-08;
     shape[155]=5.47984e-09;
     shape[156]=2.26261e-09;
     shape[157]=8.92941e-10;
     shape[158]=3.33793e-10;
     shape[159]=1.15016e-10;
     shape[160]=3.44357e-11;
     shape[161]=7.25987e-12;
     shape[162]=0;
     shape[163]=0;
     shape[164]=0;
     shape[165]=0;
     shape[166]=0;
     shape[167]=0;
     shape[168]=0;
     shape[169]=0;
     shape[170]=0;
     shape[171]=0;
     shape[172]=0;
     shape[173]=0;
     shape[174]=0;
     shape[175]=0;
     shape[176]=0;
     shape[177]=0;
     shape[178]=0;
     shape[179]=0;
     shape[180]=0;
     shape[181]=0;
     shape[182]=0;
     shape[183]=0;
     shape[184]=0;
     shape[185]=0;
     shape[186]=0;
     shape[187]=0;
     shape[188]=0;
     shape[189]=0;
     shape[190]=0;
     shape[191]=0;
     shape[192]=0;
     shape[193]=0;
     shape[194]=0;
     shape[195]=0;
     shape[196]=0;
     shape[197]=0;
     shape[198]=0;
     shape[199]=0;
     shape[200]=0;
     shape[201]=0;
     shape[202]=0;
     shape[203]=0;
     shape[204]=0;
     shape[205]=0;
     shape[206]=0;
     shape[207]=0;
     shape[208]=0;
     shape[209]=0;
     shape[210]=0;
     shape[211]=0;
     shape[212]=0;
     shape[213]=0;
     shape[214]=0;
     shape[215]=0;
     shape[216]=0;
     shape[217]=0;
     shape[218]=0;
     shape[219]=0;
     shape[220]=0;
     shape[221]=0;
     shape[222]=0;
     shape[223]=0;
     shape[224]=0;
     shape[225]=0;
     shape[226]=0;
     shape[227]=0;
     shape[228]=0;
     shape[229]=0;
     shape[230]=0;
     shape[231]=0;
     shape[232]=0;
     shape[233]=0;
     shape[234]=0;
     shape[235]=0;
     shape[236]=0;
     shape[237]=0;
     shape[238]=0;
     shape[239]=0;
     shape[240]=0;
     shape[241]=0;
     shape[242]=0;
     shape[243]=0;
     shape[244]=0;
     shape[245]=0;
     shape[246]=0;
     shape[247]=0;
     shape[248]=0;
     shape[249]=0;
     shape[250]=0;
     shape[251]=0;
     shape[252]=0;
     shape[253]=0;
     shape[254]=0;
     shape[255]=0;
     
     std::vector<float> ntmpbin(10,0.0);  // zeroing output binned pulse shape
     std::vector<float> ntmpshift(nbin,0.0);  // zeroing output shifted pulse shape  
     
     int xx = (int)x[0];
     double w1 = pars[0];
     double w2 = pars[1];
     double w3 = pars[2];
        
     TH1F *h1=new TH1F("h1","test",256,0,256);
     
     for(int i=0;i<nbin;i++){
        h1->SetBinContent(i+1,shape[i]);
     }
     
     for(int i=0;i<nbin;i++){
        if(i<w1+98.5) {
     	   ntmpshift[i]=0;
        } else {
     	   ntmpshift[i]=h1->Interpolate(i-98.5-w1);
        }
     }
     
     for(int i=0;i<nbin;i++) {
        if(i<250) {
     	   if(i<25) {
     	      ntmpbin[0]+=ntmpshift[i];
     	   } else if(i>=25&&i<50) {
     	      ntmpbin[1]+=ntmpshift[i];
     	   } else if(i>=50&&i<75) {
     	      ntmpbin[2]+=ntmpshift[i];
     	   } else if(i>=75&&i<100) {
     	      ntmpbin[3]+=ntmpshift[i];
     	   } else if(i>=100&&i<125) {
     	      ntmpbin[4]+=ntmpshift[i];
     	   } else if(i>=125&&i<150) {
     	      ntmpbin[5]+=ntmpshift[i];
     	   } else if(i>=150&&i<175) {
     	      ntmpbin[6]+=ntmpshift[i];
     	   } else if(i>=175&&i<200) {
     	      ntmpbin[7]+=ntmpshift[i];
     	   } else if(i>=200&&i<225) {
     	      ntmpbin[8]+=ntmpshift[i];
     	   } else if(i>=225&&i<250) {
     	      ntmpbin[9]+=ntmpshift[i];
     	   }
     	}
     }
     
     if( h1 ) delete h1;
     return w2*ntmpbin[xx]+w3;
  }

  double func_DoublePulse_HPDShapeData(double* x, double* pars) {
  // pulse shape componnts over a range of time 0 ns to 255 ns in 1 ns steps
     int nbin = 256;
     std::vector<float> shape(nbin,0.0);
     
     shape[0]=0;
     shape[1]=0.00033394;
     shape[2]=0.00136382;
     shape[3]=0.00345003;
     shape[4]=0.00692451;
     shape[5]=0.0114013;
     shape[6]=0.0164575;
     shape[7]=0.0216718;
     shape[8]=0.0266611;
     shape[9]=0.0311094;
     shape[10]=0.0347871;
     shape[11]=0.0375579;
     shape[12]=0.0393753;
     shape[13]=0.0402696;
     shape[14]=0.0403301;
     shape[15]=0.0396845;
     shape[16]=0.0384793;
     shape[17]=0.0368628;
     shape[18]=0.0349727;
     shape[19]=0.032928;
     shape[20]=0.0308253;
     shape[21]=0.0287377;
     shape[22]=0.0267174;
     shape[23]=0.024798;
     shape[24]=0.022999;
     shape[25]=0.0213288;
     shape[26]=0.0197884;
     shape[27]=0.0183738;
     shape[28]=0.0170779;
     shape[29]=0.0158922;
     shape[30]=0.0148077;
     shape[31]=0.0138152;
     shape[32]=0.0129064;
     shape[33]=0.012073;
     shape[34]=0.011308;
     shape[35]=0.0106047;
     shape[36]=0.00995714;
     shape[37]=0.00936001;
     shape[38]=0.00880855;
     shape[39]=0.00829848;
     shape[40]=0.007826;
     shape[41]=0.00738768;
     shape[42]=0.00698044;
     shape[43]=0.00660155;
     shape[44]=0.00624853;
     shape[45]=0.00591916;
     shape[46]=0.00561146;
     shape[47]=0.00532362;
     shape[48]=0.00505403;
     shape[49]=0.00480122;
     shape[50]=0.00456388;
     shape[51]=0.00434079;
     shape[52]=0.0041309;
     shape[53]=0.0039332;
     shape[54]=0.00374682;
     shape[55]=0.00357093;
     shape[56]=0.00340479;
     shape[57]=0.00324773;
     shape[58]=0.00309914;
     shape[59]=0.00295845;
     shape[60]=0.00282513;
     shape[61]=0.00269872;
     shape[62]=0.00257878;
     shape[63]=0.0024649;
     shape[64]=0.00235672;
     shape[65]=0.00225389;
     shape[66]=0.0021561;
     shape[67]=0.00206304;
     shape[68]=0.00197446;
     shape[69]=0.00189009;
     shape[70]=0.00180971;
     shape[71]=0.00173309;
     shape[72]=0.00166003;
     shape[73]=0.00159034;
     shape[74]=0.00152384;
     shape[75]=0.00146036;
     shape[76]=0.00139975;
     shape[77]=0.00134186;
     shape[78]=0.00128656;
     shape[79]=0.00123371;
     shape[80]=0.0011832;
     shape[81]=0.0011349;
     shape[82]=0.00108872;
     shape[83]=0.00104454;
     shape[84]=0.00100228;
     shape[85]=0.000961835;
     shape[86]=0.00092313;
     shape[87]=0.000886081;
     shape[88]=0.000850609;
     shape[89]=0.000816643;
     shape[90]=0.000784113;
     shape[91]=0.000752953;
     shape[92]=0.000723102;
     shape[93]=0.0006945;
     shape[94]=0.000667091;
     shape[95]=0.000640823;
     shape[96]=0.000615643;
     shape[97]=0.000591505;
     shape[98]=0.000568362;
     shape[99]=0.00054617;
     shape[100]=0.000524888;
     shape[101]=0.000504476;
     shape[102]=0.000484897;
     shape[103]=0.000466114;
     shape[104]=0.000448094;
     shape[105]=0.000430803;
     shape[106]=0.000414211;
     shape[107]=0.000398286;
     shape[108]=0.000383002;
     shape[109]=0.000368331;
     shape[110]=0.000354248;
     shape[111]=0.000340726;
     shape[112]=0.000327743;
     shape[113]=0.000315276;
     shape[114]=0.000303304;
     shape[115]=0.000291806;
     shape[116]=0.000280762;
     shape[117]=0.000270153;
     shape[118]=0.000259962;
     shape[119]=0.000250171;
     shape[120]=0.000240764;
     shape[121]=0.000231725;
     shape[122]=0.000223038;
     shape[123]=0.000214691;
     shape[124]=0.000206667;
     shape[125]=0.000198956;
     shape[126]=0.000191543;
     shape[127]=0.000184416;
     shape[128]=0.000177565;
     shape[129]=0.000170978;
     shape[130]=0.000164644;
     shape[131]=0.000158554;
     shape[132]=0.000152697;
     shape[133]=0.000147064;
     shape[134]=0.000141646;
     shape[135]=0.000136435;
     shape[136]=0.000131423;
     shape[137]=0.000126601;
     shape[138]=0.000121962;
     shape[139]=0.000117498;
     shape[140]=0.000113204;
     shape[141]=0.000109071;
     shape[142]=0.000105095;
     shape[143]=0.000101268;
     shape[144]=9.75848e-05;
     shape[145]=9.404e-05;
     shape[146]=9.0628e-05;
     shape[147]=8.73436e-05;
     shape[148]=8.4182e-05;
     shape[149]=8.11383e-05;
     shape[150]=7.8208e-05;
     shape[151]=7.53866e-05;
     shape[152]=7.26701e-05;
     shape[153]=7.00543e-05;
     shape[154]=6.75353e-05;
     shape[155]=6.51096e-05;
     shape[156]=6.27734e-05;
     shape[157]=6.05233e-05;
     shape[158]=5.83562e-05;
     shape[159]=5.62687e-05;
     shape[160]=5.42579e-05;
     shape[161]=5.23209e-05;
     shape[162]=5.04549e-05;
     shape[163]=4.86571e-05;
     shape[164]=4.69251e-05;
     shape[165]=4.52562e-05;
     shape[166]=4.36482e-05;
     shape[167]=4.20987e-05;
     shape[168]=4.06056e-05;
     shape[169]=3.91667e-05;
     shape[170]=3.778e-05;
     shape[171]=3.64436e-05;
     shape[172]=3.51555e-05;
     shape[173]=3.3914e-05;
     shape[174]=3.27173e-05;
     shape[175]=3.14624e-05;
     shape[176]=3.00317e-05;
     shape[177]=2.83024e-05;
     shape[178]=2.61543e-05;
     shape[179]=2.36806e-05;
     shape[180]=2.09939e-05;
     shape[181]=1.82145e-05;
     shape[182]=1.5459e-05;
     shape[183]=1.28302e-05;
     shape[184]=1.041e-05;
     shape[185]=8.25538e-06;
     shape[186]=6.3974e-06;
     shape[187]=4.84374e-06;
     shape[188]=3.58271e-06;
     shape[189]=2.58848e-06;
     shape[190]=1.82658e-06;
     shape[191]=1.25879e-06;
     shape[192]=8.47151e-07;
     shape[193]=5.56713e-07;
     shape[194]=3.57223e-07;
     shape[195]=2.23802e-07;
     shape[196]=1.36893e-07;
     shape[197]=8.17475e-08;
     shape[198]=4.76566e-08;
     shape[199]=2.71208e-08;
     shape[200]=1.50656e-08;
     shape[201]=8.16834e-09;
     shape[202]=4.32189e-09;
     shape[203]=2.2309e-09;
     shape[204]=1.12281e-09;
     shape[205]=5.50335e-10;
     shape[206]=2.62018e-10;
     shape[207]=1.20477e-10;
     shape[208]=5.27642e-11;
     shape[209]=2.12197e-11;
     shape[210]=7.34902e-12;
     shape[211]=1.76216e-12;
     shape[212]=0;
     shape[213]=0;
     shape[214]=0;
     shape[215]=0;
     shape[216]=0;
     shape[217]=0;
     shape[218]=0;
     shape[219]=0;
     shape[220]=0;
     shape[221]=0;
     shape[222]=0;
     shape[223]=0;
     shape[224]=0;
     shape[225]=0;
     shape[226]=0;
     shape[227]=0;
     shape[228]=0;
     shape[229]=0;
     shape[230]=0;
     shape[231]=0;
     shape[232]=0;
     shape[233]=0;
     shape[234]=0;
     shape[235]=0;
     shape[236]=0;
     shape[237]=0;
     shape[238]=0;
     shape[239]=0;
     shape[240]=0;
     shape[241]=0;
     shape[242]=0;
     shape[243]=0;
     shape[244]=0;
     shape[245]=0;
     shape[246]=0;
     shape[247]=0;
     shape[248]=0;
     shape[249]=0;
     shape[250]=0;
     shape[251]=0;
     shape[252]=0;
     shape[253]=0;
     shape[254]=0;
     shape[255]=0;
     
     std::vector<float> ntmpbin(10,0.0);  // zeroing output binned pulse shape
     std::vector<float> ntmpbin2(10,0.0);  // zeroing output binned pulse shape
     std::vector<float> ntmpshift(nbin,0.0);  // zeroing output shifted pulse shape  
     std::vector<float> ntmpshift2(nbin,0.0);  // zeroing output shifted pulse shape
     
     int xx = (int)x[0];
     double w1 = pars[0];
     double w2 = pars[1];
     double w3 = pars[2];
     double w4 = pars[3];
     double w5 = pars[4];
        
     TH1F *h1=new TH1F("h1","test",256,0,256);
     
     for(int i=0;i<nbin;i++){
        h1->SetBinContent(i+1,shape[i]);
     }
     
     for(int i=0;i<nbin;i++){
        if(i<w1+98.5) {
     	   ntmpshift[i]=0;
        } else {
     	   ntmpshift[i]=h1->Interpolate(i-98.5-w1);
        }
     }
     
     for(int i=0;i<nbin;i++){
        if(i<w2+98.5) {
     	   ntmpshift2[i]=0;
        } else {
     	   ntmpshift2[i]=h1->Interpolate(i-98.5-w2);
        }
     }
     
     for(int i=0;i<nbin;i++) {
        if(i<250) {
     	   if(i<25) {
     	      ntmpbin[0]+=ntmpshift[i];
     	      ntmpbin2[0]+=ntmpshift2[i];
     	   } else if(i>=25&&i<50) {
     	      ntmpbin[1]+=ntmpshift[i];
     	      ntmpbin2[1]+=ntmpshift2[i];
     	   } else if(i>=50&&i<75) {
     	      ntmpbin[2]+=ntmpshift[i];
     	      ntmpbin2[2]+=ntmpshift2[i];
     	   } else if(i>=75&&i<100) {
     	      ntmpbin[3]+=ntmpshift[i];
     	      ntmpbin2[3]+=ntmpshift2[i];
     	   } else if(i>=100&&i<125) {
     	      ntmpbin[4]+=ntmpshift[i];
     	      ntmpbin2[4]+=ntmpshift2[i];
     	   } else if(i>=125&&i<150) {
     	      ntmpbin[5]+=ntmpshift[i];
     	      ntmpbin2[5]+=ntmpshift2[i];
     	   } else if(i>=150&&i<175) {
     	      ntmpbin[6]+=ntmpshift[i];
     	      ntmpbin2[6]+=ntmpshift2[i];
     	   } else if(i>=175&&i<200) {
     	      ntmpbin[7]+=ntmpshift[i];
     	      ntmpbin2[7]+=ntmpshift2[i];
     	   } else if(i>=200&&i<225) {
     	      ntmpbin[8]+=ntmpshift[i];
     	      ntmpbin2[8]+=ntmpshift2[i];
     	   } else if(i>=225&&i<250) {
     	      ntmpbin[9]+=ntmpshift[i];
     	      ntmpbin2[9]+=ntmpshift2[i];
     	   }
     	}
     }
     
     if( h1 ) delete h1;
     return w3*ntmpbin[xx]+w4*ntmpbin2[xx]+w5;
  }

  double func_DoublePulse_HPDShapeMC(double* x, double* pars) {
  // pulse shape componnts over a range of time 0 ns to 255 ns in 1 ns steps
     int nbin = 256;
     std::vector<float> shape(nbin,0.0);
     
     shape[0]=0;
     shape[1]=0.000498985;
     shape[2]=0.00201622;
     shape[3]=0.00502753;
     shape[4]=0.00991548;
     shape[5]=0.0159424;
     shape[6]=0.0223257;
     shape[7]=0.0283523;
     shape[8]=0.0334672;
     shape[9]=0.0373233;
     shape[10]=0.0397869;
     shape[11]=0.0409066;
     shape[12]=0.0408609;
     shape[13]=0.0398981;
     shape[14]=0.0382843;
     shape[15]=0.0362642;
     shape[16]=0.0340391;
     shape[17]=0.0317595;
     shape[18]=0.0295275;
     shape[19]=0.0274053;
     shape[20]=0.0254252;
     shape[21]=0.0235992;
     shape[22]=0.0219269;
     shape[23]=0.0204009;
     shape[24]=0.0190103;
     shape[25]=0.0177433;
     shape[26]=0.0165881;
     shape[27]=0.0155337;
     shape[28]=0.0145699;
     shape[29]=0.0136876;
     shape[30]=0.0128785;
     shape[31]=0.0121354;
     shape[32]=0.0114516;
     shape[33]=0.0108214;
     shape[34]=0.0102396;
     shape[35]=0.00970139;
     shape[36]=0.00920274;
     shape[37]=0.00873989;
     shape[38]=0.00830953;
     shape[39]=0.00790869;
     shape[40]=0.00753469;
     shape[41]=0.00718516;
     shape[42]=0.00685795;
     shape[43]=0.00655115;
     shape[44]=0.00626302;
     shape[45]=0.00599202;
     shape[46]=0.00573674;
     shape[47]=0.00549594;
     shape[48]=0.00526847;
     shape[49]=0.00505331;
     shape[50]=0.00484954;
     shape[51]=0.00465631;
     shape[52]=0.00447286;
     shape[53]=0.00429852;
     shape[54]=0.00413264;
     shape[55]=0.00397466;
     shape[56]=0.00382407;
     shape[57]=0.00368038;
     shape[58]=0.00354316;
     shape[59]=0.00341202;
     shape[60]=0.00328659;
     shape[61]=0.00316654;
     shape[62]=0.00305157;
     shape[63]=0.00294138;
     shape[64]=0.00283572;
     shape[65]=0.00273435;
     shape[66]=0.00263704;
     shape[67]=0.00254359;
     shape[68]=0.00245379;
     shape[69]=0.00236749;
     shape[70]=0.0022845;
     shape[71]=0.00220466;
     shape[72]=0.00212785;
     shape[73]=0.0020539;
     shape[74]=0.00198271;
     shape[75]=0.00191415;
     shape[76]=0.0018481;
     shape[77]=0.00178445;
     shape[78]=0.00172312;
     shape[79]=0.00166399;
     shape[80]=0.00160699;
     shape[81]=0.00155202;
     shape[82]=0.001499;
     shape[83]=0.00144786;
     shape[84]=0.00139852;
     shape[85]=0.00135092;
     shape[86]=0.00130498;
     shape[87]=0.00126065;
     shape[88]=0.00121787;
     shape[89]=0.00117656;
     shape[90]=0.00113669;
     shape[91]=0.0010982;
     shape[92]=0.00106104;
     shape[93]=0.00102515;
     shape[94]=0.000990499;
     shape[95]=0.000957035;
     shape[96]=0.000924718;
     shape[97]=0.000893506;
     shape[98]=0.00086336;
     shape[99]=0.000834243;
     shape[100]=0.000806118;
     shape[101]=0.00077895;
     shape[102]=0.000752705;
     shape[103]=0.000727352;
     shape[104]=0.00070286;
     shape[105]=0.000679198;
     shape[106]=0.000656338;
     shape[107]=0.000634253;
     shape[108]=0.000612914;
     shape[109]=0.000592298;
     shape[110]=0.000572378;
     shape[111]=0.000553131;
     shape[112]=0.000534535;
     shape[113]=0.000516565;
     shape[114]=0.000499203;
     shape[115]=0.000482425;
     shape[116]=0.000466214;
     shape[117]=0.000450549;
     shape[118]=0.000435411;
     shape[119]=0.000420784;
     shape[120]=0.000406649;
     shape[121]=0.00039299;
     shape[122]=0.000379791;
     shape[123]=0.000367036;
     shape[124]=0.00035471;
     shape[125]=0.000342798;
     shape[126]=0.000331288;
     shape[127]=0.000320164;
     shape[128]=0.000309414;
     shape[129]=0.000299026;
     shape[130]=0.000288987;
     shape[131]=0.000279285;
     shape[132]=0.00026991;
     shape[133]=0.000260849;
     shape[134]=0.000252093;
     shape[135]=0.000243631;
     shape[136]=0.000235453;
     shape[137]=0.00022755;
     shape[138]=0.000219912;
     shape[139]=0.000212531;
     shape[140]=0.000205398;
     shape[141]=0.000198504;
     shape[142]=0.000191842;
     shape[143]=0.000185404;
     shape[144]=0.000179181;
     shape[145]=0.000173168;
     shape[146]=0.000167356;
     shape[147]=0.00016174;
     shape[148]=0.000156312;
     shape[149]=0.000151066;
     shape[150]=0.000145997;
     shape[151]=0.000141098;
     shape[152]=0.000136363;
     shape[153]=0.000131787;
     shape[154]=0.000127364;
     shape[155]=0.00012309;
     shape[156]=0.00011896;
     shape[157]=0.000114968;
     shape[158]=0.00011111;
     shape[159]=0.000107381;
     shape[160]=0.000103778;
     shape[161]=0.000100296;
     shape[162]=9.69302e-05;
     shape[163]=9.36776e-05;
     shape[164]=9.05342e-05;
     shape[165]=8.74963e-05;
     shape[166]=8.45603e-05;
     shape[167]=8.17229e-05;
     shape[168]=7.89806e-05;
     shape[169]=7.63304e-05;
     shape[170]=7.37691e-05;
     shape[171]=7.12938e-05;
     shape[172]=6.89016e-05;
     shape[173]=6.65896e-05;
     shape[174]=6.43552e-05;
     shape[175]=6.18403e-05;
     shape[176]=5.86522e-05;
     shape[177]=5.44093e-05;
     shape[178]=4.87793e-05;
     shape[179]=4.22159e-05;
     shape[180]=3.52384e-05;
     shape[181]=2.83509e-05;
     shape[182]=2.19736e-05;
     shape[183]=1.63998e-05;
     shape[184]=1.17825e-05;  
     shape[185]=8.14666e-06;
     shape[186]=5.41959e-06;
     shape[187]=3.46829e-06;
     shape[188]=2.13479e-06;
     shape[189]=1.26364e-06;
     shape[190]=7.19224e-07;
     shape[191]=3.93576e-07;
     shape[192]=2.07047e-07;
     shape[193]=1.04697e-07;
     shape[194]=5.08839e-08;
     shape[195]=2.37647e-08;
     shape[196]=1.0663e-08;
     shape[197]=4.59391e-09;
     shape[198]=1.89802e-09;
     shape[199]=7.49603e-10;
     shape[200]=2.80472e-10;
     shape[201]=9.67625e-11;
     shape[202]=2.90206e-11;
     shape[203]=6.1327e-12;
     shape[204]=0;
     shape[205]=0;
     shape[206]=0;
     shape[207]=0;
     shape[208]=0;
     shape[209]=0;
     shape[210]=0;
     shape[211]=0;
     shape[212]=0;
     shape[213]=0;
     shape[214]=0;
     shape[215]=0;
     shape[216]=0;
     shape[217]=0;
     shape[218]=0;
     shape[219]=0;
     shape[220]=0;
     shape[221]=0;
     shape[222]=0;
     shape[223]=0;
     shape[224]=0;
     shape[225]=0;
     shape[226]=0;
     shape[227]=0;
     shape[228]=0;
     shape[229]=0;
     shape[230]=0;
     shape[231]=0;
     shape[232]=0;
     shape[233]=0;
     shape[234]=0;
     shape[235]=0;
     shape[236]=0;
     shape[237]=0;
     shape[238]=0;
     shape[239]=0;
     shape[240]=0;
     shape[241]=0;
     shape[242]=0;
     shape[243]=0;
     shape[244]=0;
     shape[245]=0;
     shape[246]=0;
     shape[247]=0;
     shape[248]=0;
     shape[249]=0;
     shape[250]=0;
     shape[251]=0;
     shape[252]=0;
     shape[253]=0;
     shape[254]=0;
     shape[255]=0;
     
     std::vector<float> ntmpbin(10,0.0);  // zeroing output binned pulse shape
     std::vector<float> ntmpbin2(10,0.0);  // zeroing output binned pulse shape
     std::vector<float> ntmpshift(nbin,0.0);  // zeroing output shifted pulse shape  
     std::vector<float> ntmpshift2(nbin,0.0);  // zeroing output shifted pulse shape
     
     int xx = (int)x[0];
     double w1 = pars[0];
     double w2 = pars[1];
     double w3 = pars[2];
     double w4 = pars[3];
     double w5 = pars[4];
        
     TH1F *h1=new TH1F("h1","test",256,0,256);
     
     for(int i=0;i<nbin;i++){
        h1->SetBinContent(i+1,shape[i]);
     }
     
     for(int i=0;i<nbin;i++){
        if(i<w1+98.5) {
       	   ntmpshift[i]=0;
        } else {
     	   ntmpshift[i]=h1->Interpolate(i-98.5-w1);
        }
     }
     
     for(int i=0;i<nbin;i++){
        if(i<w2+98.5) {
     	   ntmpshift2[i]=0;
        } else {
     	   ntmpshift2[i]=h1->Interpolate(i-98.5-w2);
        }
    }
     
    for(int i=0;i<nbin;i++) {
       if(i<250) {
     	  if(i<25) {
     	     ntmpbin[0]+=ntmpshift[i];
     	     ntmpbin2[0]+=ntmpshift2[i];
     	  } else if(i>=25&&i<50) {
     	     ntmpbin[1]+=ntmpshift[i];
     	     ntmpbin2[1]+=ntmpshift2[i];
     	  } else if(i>=50&&i<75) {
     	     ntmpbin[2]+=ntmpshift[i];
     	     ntmpbin2[2]+=ntmpshift2[i];
     	  } else if(i>=75&&i<100) {
     	     ntmpbin[3]+=ntmpshift[i];
     	     ntmpbin2[3]+=ntmpshift2[i];
     	  } else if(i>=100&&i<125) {
     	     ntmpbin[4]+=ntmpshift[i];
     	     ntmpbin2[4]+=ntmpshift2[i];
     	  } else if(i>=125&&i<150) {
     	     ntmpbin[5]+=ntmpshift[i];
     	     ntmpbin2[5]+=ntmpshift2[i];
     	  } else if(i>=150&&i<175) {
     	     ntmpbin[6]+=ntmpshift[i];
     	     ntmpbin2[6]+=ntmpshift2[i];
     	  } else if(i>=175&&i<200) {
     	     ntmpbin[7]+=ntmpshift[i];
     	     ntmpbin2[7]+=ntmpshift2[i];
     	  } else if(i>=200&&i<225) {
     	     ntmpbin[8]+=ntmpshift[i];
     	     ntmpbin2[8]+=ntmpshift2[i];
     	  } else if(i>=225&&i<250) {
     	     ntmpbin[9]+=ntmpshift[i];
     	     ntmpbin2[9]+=ntmpshift2[i];
     	  }
       }
    }
   
    if( h1 ) delete h1;  
    return w3*ntmpbin[xx]+w4*ntmpbin2[xx]+w5;
  }
   
}


HBHERecHit HcalSimpleRecAlgo::reconstruct(const HBHEDataFrame& digi, int first, int toadd, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  return HcalSimpleRecAlgoImpl::reco<HBHEDataFrame,HBHERecHit>(digi,coder,calibs,
							       first,toadd,correctForTimeslew_, correctForPulse_,
							       pulseCorr_->get(digi.id(), toadd, phaseNS_),
							       HcalTimeSlew::Medium,
                                                               runnum_, setLeakCorrection_,
                                                               hbhePileupCorr_.get(),
                                                               bunchCrossingInfo_, lenBunchCrossingInfo_, puCorrMethod_);
}


HORecHit HcalSimpleRecAlgo::reconstruct(const HODataFrame& digi, int first, int toadd, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  return HcalSimpleRecAlgoImpl::reco<HODataFrame,HORecHit>(digi,coder,calibs,
							   first,toadd,correctForTimeslew_,correctForPulse_,
							   pulseCorr_->get(digi.id(), toadd, phaseNS_),
							   HcalTimeSlew::Slow,
                                                           runnum_, false, hoPileupCorr_.get(),
                                                           bunchCrossingInfo_, lenBunchCrossingInfo_, puCorrMethod_);
}


HcalCalibRecHit HcalSimpleRecAlgo::reconstruct(const HcalCalibDataFrame& digi, int first, int toadd, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  return HcalSimpleRecAlgoImpl::reco<HcalCalibDataFrame,HcalCalibRecHit>(digi,coder,calibs,
									 first,toadd,correctForTimeslew_,correctForPulse_,
									 pulseCorr_->get(digi.id(), toadd, phaseNS_),
									 HcalTimeSlew::Fast,
                                                                         runnum_, false, 0,
                                                                         bunchCrossingInfo_, lenBunchCrossingInfo_, puCorrMethod_);
}


HBHERecHit HcalSimpleRecAlgo::reconstructHBHEUpgrade(const HcalUpgradeDataFrame& digi, int first, int toadd, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  HBHERecHit result = HcalSimpleRecAlgoImpl::reco<HcalUpgradeDataFrame,HBHERecHit>(digi, coder, calibs,
                                                                                   first, toadd, correctForTimeslew_, correctForPulse_,
                                                                                   pulseCorr_->get(digi.id(), toadd, phaseNS_),
                                                                                   HcalTimeSlew::Medium, 0, false,
                                                                                   hbhePileupCorr_.get(),
                                                                                   bunchCrossingInfo_, lenBunchCrossingInfo_, puCorrMethod_);
  HcalTDCReco tdcReco;
  tdcReco.reconstruct(digi, result);
  return result;
}


HFRecHit HcalSimpleRecAlgo::reconstruct(const HFDataFrame& digi,
                                        const int first,
                                        const int toadd,
                                        const HcalCoder& coder,
                                        const HcalCalibrations& calibs) const
{
  const HcalPulseContainmentCorrection* corr = pulseCorr_->get(digi.id(), toadd, phaseNS_);

  double amp_fC, ampl, uncorr_ampl, maxA;
  int nRead, maxI;
  bool leakCorrApplied;
  float t0, t2;

  HcalSimpleRecAlgoImpl::removePileup(digi, coder, calibs, first, toadd,
                                      correctForPulse_, corr, hfPileupCorr_.get(),
                                      bunchCrossingInfo_, lenBunchCrossingInfo_,
                                      &maxA, &ampl, &uncorr_ampl, &amp_fC, &nRead,
                                      &maxI, &leakCorrApplied, &t0, &t2);

  float time=-9999.f;
  if (maxI > 0 && maxI < (nRead - 1))
      time = HcalSimpleRecAlgoImpl::recoHFTime(digi,maxI,amp_fC,correctForTimeslew_,maxA,t0,t2) -
             calibs.timecorr();

  HFRecHit rh(digi.id(),ampl,time);
  setRawEnergy(rh, static_cast<float>(uncorr_ampl));
  return rh;
}


// NB: Upgrade HFRecHit method content is just the same as regular  HFRecHit
//     with one exclusion: double time (second is dummy) in constructor 
HFRecHit HcalSimpleRecAlgo::reconstructHFUpgrade(const HcalUpgradeDataFrame& digi,
                                                 const int first,
                                                 const int toadd,
                                                 const HcalCoder& coder,
                                                 const HcalCalibrations& calibs) const
{
  const HcalPulseContainmentCorrection* corr = pulseCorr_->get(digi.id(), toadd, phaseNS_);

  double amp_fC, ampl, uncorr_ampl, maxA;
  int nRead, maxI;
  bool leakCorrApplied;
  float t0, t2;

  HcalSimpleRecAlgoImpl::removePileup(digi, coder, calibs, first, toadd,
                                      correctForPulse_, corr, hfPileupCorr_.get(),
                                      bunchCrossingInfo_, lenBunchCrossingInfo_,
                                      &maxA, &ampl, &uncorr_ampl, &amp_fC, &nRead,
                                      &maxI, &leakCorrApplied, &t0, &t2);

  float time=-9999.f;
  if (maxI > 0 && maxI < (nRead - 1))
      time = HcalSimpleRecAlgoImpl::recoHFTime(digi,maxI,amp_fC,correctForTimeslew_,maxA,t0,t2) -
             calibs.timecorr();

  HFRecHit rh(digi.id(),ampl,time); // new RecHit gets second time = 0.
  setRawEnergy(rh, static_cast<float>(uncorr_ampl));
  return rh;
}


/// Ugly hack to apply energy corrections to some HB- cells
float eCorr(int ieta, int iphi, double energy, int runnum) {
// return energy correction factor for HBM channels 
// iphi=6 ieta=(-1,-15) and iphi=32 ieta=(-1,-7)
// I.Vodopianov 28 Feb. 2011
  static const float low32[7]  = {0.741,0.721,0.730,0.698,0.708,0.751,0.861};
  static const float high32[7] = {0.973,0.925,0.900,0.897,0.950,0.935,1};
  static const float low6[15]  = {0.635,0.623,0.670,0.633,0.644,0.648,0.600,
				  0.570,0.595,0.554,0.505,0.513,0.515,0.561,0.579};
  static const float high6[15] = {0.875,0.937,0.942,0.900,0.922,0.925,0.901,
				  0.850,0.852,0.818,0.731,0.717,0.782,0.853,0.778};

  
  double slope, mid, en;
  double corr = 1.0;

  if (!(iphi==6 && ieta<0 && ieta>-16) && !(iphi==32 && ieta<0 && ieta>-8)) 
    return corr;

  int jeta = -ieta-1;
  double xeta = (double) ieta;
  if (energy > 0.) en=energy;
  else en = 0.;

  if (iphi == 32) {
    slope = 0.2272;
    mid = 17.14 + 0.7147*xeta;
    if (en > 100.) corr = high32[jeta];
    else corr = low32[jeta]+(high32[jeta]-low32[jeta])/(1.0+exp(-(en-mid)*slope));
  }
  else if (iphi == 6 && runnum < 216091 ) {
    slope = 0.1956;
    mid = 15.96 + 0.3075*xeta;
    if (en > 100.0) corr = high6[jeta];
    else corr = low6[jeta]+(high6[jeta]-low6[jeta])/(1.0+exp(-(en-mid)*slope));
  }

  //  std::cout << "HBHE cell:  ieta, iphi = " << ieta << "  " << iphi 
  //	    << "  ->  energy = " << en << "   corr = " << corr << std::endl;

  return corr;
}


// Actual leakage (to pre-sample) correction 
float leakCorr(double energy) {
  double corr = 1.0;
  return corr;
}


// timeshift implementation

static const float wpksamp0_hbheho = 0.5;
static const int   num_bins_hbheho = 61;

static const float actual_ns_hbheho[num_bins_hbheho] = {
-5.44000, // 0.500, 0.000-0.017
-4.84250, // 0.517, 0.017-0.033
-4.26500, // 0.533, 0.033-0.050
-3.71000, // 0.550, 0.050-0.067
-3.18000, // 0.567, 0.067-0.083
-2.66250, // 0.583, 0.083-0.100
-2.17250, // 0.600, 0.100-0.117
-1.69000, // 0.617, 0.117-0.133
-1.23000, // 0.633, 0.133-0.150
-0.78000, // 0.650, 0.150-0.167
-0.34250, // 0.667, 0.167-0.183
 0.08250, // 0.683, 0.183-0.200
 0.50250, // 0.700, 0.200-0.217
 0.90500, // 0.717, 0.217-0.233
 1.30500, // 0.733, 0.233-0.250
 1.69500, // 0.750, 0.250-0.267
 2.07750, // 0.767, 0.267-0.283
 2.45750, // 0.783, 0.283-0.300
 2.82500, // 0.800, 0.300-0.317
 3.19250, // 0.817, 0.317-0.333
 3.55750, // 0.833, 0.333-0.350
 3.91750, // 0.850, 0.350-0.367
 4.27500, // 0.867, 0.367-0.383
 4.63000, // 0.883, 0.383-0.400
 4.98500, // 0.900, 0.400-0.417
 5.33750, // 0.917, 0.417-0.433
 5.69500, // 0.933, 0.433-0.450
 6.05000, // 0.950, 0.450-0.467
 6.40500, // 0.967, 0.467-0.483
 6.77000, // 0.983, 0.483-0.500
 7.13500, // 1.000, 0.500-0.517
 7.50000, // 1.017, 0.517-0.533
 7.88250, // 1.033, 0.533-0.550
 8.26500, // 1.050, 0.550-0.567
 8.66000, // 1.067, 0.567-0.583
 9.07000, // 1.083, 0.583-0.600
 9.48250, // 1.100, 0.600-0.617
 9.92750, // 1.117, 0.617-0.633
10.37750, // 1.133, 0.633-0.650
10.87500, // 1.150, 0.650-0.667
11.38000, // 1.167, 0.667-0.683
11.95250, // 1.183, 0.683-0.700
12.55000, // 1.200, 0.700-0.717
13.22750, // 1.217, 0.717-0.733
13.98500, // 1.233, 0.733-0.750
14.81500, // 1.250, 0.750-0.767
15.71500, // 1.267, 0.767-0.783
16.63750, // 1.283, 0.783-0.800
17.53750, // 1.300, 0.800-0.817
18.38500, // 1.317, 0.817-0.833
19.16500, // 1.333, 0.833-0.850
19.89750, // 1.350, 0.850-0.867
20.59250, // 1.367, 0.867-0.883
21.24250, // 1.383, 0.883-0.900
21.85250, // 1.400, 0.900-0.917
22.44500, // 1.417, 0.917-0.933
22.99500, // 1.433, 0.933-0.950
23.53250, // 1.450, 0.950-0.967
24.03750, // 1.467, 0.967-0.983
24.53250, // 1.483, 0.983-1.000
25.00000  // 1.500, 1.000-1.017 - keep for interpolation
};

float timeshift_ns_hbheho(float wpksamp) {
  float flx = (num_bins_hbheho-1)*(wpksamp - wpksamp0_hbheho);
  int index = (int)flx;
  float yval;

  if      (index <    0)               return actual_ns_hbheho[0];
  else if (index >= num_bins_hbheho-1) return actual_ns_hbheho[num_bins_hbheho-1];

  // else interpolate:
  float y1 = actual_ns_hbheho[index];
  float y2 = actual_ns_hbheho[index+1];

  yval = y1 + (y2-y1)*(flx-(float)index);

  return yval;
}

static const int   num_bins_hf = 101;
static const float wpksamp0_hf = 0.5;

static const float actual_ns_hf[num_bins_hf] = {
 0.00250, // 0.000-0.010
 0.04500, // 0.010-0.020
 0.08750, // 0.020-0.030
 0.13000, // 0.030-0.040
 0.17250, // 0.040-0.050
 0.21500, // 0.050-0.060
 0.26000, // 0.060-0.070
 0.30250, // 0.070-0.080
 0.34500, // 0.080-0.090
 0.38750, // 0.090-0.100
 0.42750, // 0.100-0.110
 0.46000, // 0.110-0.120
 0.49250, // 0.120-0.130
 0.52500, // 0.130-0.140
 0.55750, // 0.140-0.150
 0.59000, // 0.150-0.160
 0.62250, // 0.160-0.170
 0.65500, // 0.170-0.180
 0.68750, // 0.180-0.190
 0.72000, // 0.190-0.200
 0.75250, // 0.200-0.210
 0.78500, // 0.210-0.220
 0.81750, // 0.220-0.230
 0.85000, // 0.230-0.240
 0.88250, // 0.240-0.250
 0.91500, // 0.250-0.260
 0.95500, // 0.260-0.270
 0.99250, // 0.270-0.280
 1.03250, // 0.280-0.290
 1.07000, // 0.290-0.300
 1.10750, // 0.300-0.310
 1.14750, // 0.310-0.320
 1.18500, // 0.320-0.330
 1.22500, // 0.330-0.340
 1.26250, // 0.340-0.350
 1.30000, // 0.350-0.360
 1.34000, // 0.360-0.370
 1.37750, // 0.370-0.380
 1.41750, // 0.380-0.390
 1.48750, // 0.390-0.400
 1.55750, // 0.400-0.410
 1.62750, // 0.410-0.420
 1.69750, // 0.420-0.430
 1.76750, // 0.430-0.440
 1.83750, // 0.440-0.450
 1.90750, // 0.450-0.460
 2.06750, // 0.460-0.470
 2.23250, // 0.470-0.480
 2.40000, // 0.480-0.490
 2.82250, // 0.490-0.500
 3.81000, // 0.500-0.510
 6.90500, // 0.510-0.520
 8.99250, // 0.520-0.530
10.50000, // 0.530-0.540
11.68250, // 0.540-0.550
12.66250, // 0.550-0.560
13.50250, // 0.560-0.570
14.23750, // 0.570-0.580
14.89750, // 0.580-0.590
15.49000, // 0.590-0.600
16.03250, // 0.600-0.610
16.53250, // 0.610-0.620
17.00000, // 0.620-0.630
17.44000, // 0.630-0.640
17.85250, // 0.640-0.650
18.24000, // 0.650-0.660
18.61000, // 0.660-0.670
18.96750, // 0.670-0.680
19.30500, // 0.680-0.690
19.63000, // 0.690-0.700
19.94500, // 0.700-0.710
20.24500, // 0.710-0.720
20.54000, // 0.720-0.730
20.82250, // 0.730-0.740
21.09750, // 0.740-0.750
21.37000, // 0.750-0.760
21.62750, // 0.760-0.770
21.88500, // 0.770-0.780
22.13000, // 0.780-0.790
22.37250, // 0.790-0.800
22.60250, // 0.800-0.810
22.83000, // 0.810-0.820
23.04250, // 0.820-0.830
23.24500, // 0.830-0.840
23.44250, // 0.840-0.850
23.61000, // 0.850-0.860
23.77750, // 0.860-0.870
23.93500, // 0.870-0.880
24.05500, // 0.880-0.890
24.17250, // 0.890-0.900
24.29000, // 0.900-0.910
24.40750, // 0.910-0.920
24.48250, // 0.920-0.930
24.55500, // 0.930-0.940
24.62500, // 0.940-0.950
24.69750, // 0.950-0.960
24.77000, // 0.960-0.970
24.84000, // 0.970-0.980
24.91250, // 0.980-0.990
24.95500, // 0.990-1.000
24.99750, // 1.000-1.010 - keep for interpolation
};

float timeshift_ns_hf(float wpksamp) {
  float flx = (num_bins_hf-1)*(wpksamp-wpksamp0_hf);
  int index = (int)flx;
  float yval;
  
  if      (index <  0)             return actual_ns_hf[0];
  else if (index >= num_bins_hf-1) return actual_ns_hf[num_bins_hf-1];

  // else interpolate:
  float y1       = actual_ns_hf[index];
  float y2       = actual_ns_hf[index+1];

  // float delta_x  = 1/(float)num_bins_hf;
  // yval = y1 + (y2-y1)*(flx-(float)index)/delta_x;

  yval = y1 + (y2-y1)*(flx-(float)index);
  return yval;
}

