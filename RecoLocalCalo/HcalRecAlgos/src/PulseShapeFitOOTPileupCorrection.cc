#include <climits>
#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFitOOTPileupCorrection.h"
#include "RecoLocalCalo/HcalRecAlgos/src/Asa047.h"

namespace FitterFuncs{

   int cntNANinfit;
   double psFit_x[10], psFit_y[10], psFit_erry[10]; 
  // since we know we are counting in nanoseconds
  // we don't need to do an expensive finding of the bin
  // simply take floor(x) and determine if bin center is above or below
  // bin center is just bin + 0.5, inputs bins are 1ns wide
   float fast_interpolate(double x, const std::array<float,256>& h1) {
      if( x != x ){
         cntNANinfit ++;
         return h1[255];
      }

      const int bin = (int)x;

      if( x < 0.5 ) return h1[0];
      else if ( x > 255.5 ) return h1[255];

      const int bin_0 = ( x < bin+0.5 ? bin-1 : bin );
      const int bin_1 = ( x < bin+0.5 ? bin : bin+1 );

      const float slope = (h1[bin_1] - h1[bin_0])/(bin_1-bin_0);
      return h1[bin_0] + (x-bin_0-0.5f)*slope;
   }

   std::array<float,10> funcHPDShape(const std::vector<double>& pars,
                                     const std::array<float,256>& h1_single) {
    // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
      constexpr int ns_per_bx = 25;
      constexpr int num_ns = 250;
      constexpr int num_bx = num_ns/ns_per_bx;

    // zeroing output binned pulse shape
      std::array<float,num_bx> ntmpbin{ {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f} };

      for(int i=0;i < num_ns; ++i) {
         const float offset = i - 98.5f - pars[0]; // where does 98.5 come from?
         const float shifted_pulse1 = (offset < 0.0f ? 0.0f : fast_interpolate(offset,h1_single));
         ntmpbin[i/ns_per_bx] += shifted_pulse1;
      }
    // now we use ntmpbin to record the final pulse shape
      for(int i=0; i < num_bx; ++i) {
         ntmpbin[i] = pars[1]*ntmpbin[i] + pars[2];
      }
      return ntmpbin;
   }

   std::array<float,10> func_DoublePulse_HPDShape(const std::vector<double>& pars,
                                                  const std::array<float,256>& h1_double) {
    // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
      constexpr int ns_per_bx = 25;
      constexpr int num_ns = 250;
      constexpr int num_bx = num_ns/ns_per_bx;

    // zeroing output binned pulse shape
      std::array<float,num_bx> ntmpbin{ {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f} };
      std::array<float,num_bx> ntmpbin2{ {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f} };

      for(int i=0;i < num_ns;++i) {
         const float offset1 = i - 98.5 - pars[0]; // where does 98.5 come from?
         const float offset2 = i - 98.5 - pars[1];

         ntmpbin[i/ns_per_bx] += (offset1 < 0.0f ? 0.0f : fast_interpolate(offset1,h1_double));
         ntmpbin2[i/ns_per_bx] += (offset2 < 0.0f ? 0.0f : fast_interpolate(offset2,h1_double));
      }
    // now we use ntmpbin to record the final pulse shape
      for(int i=0; i < num_bx; ++i) {
         ntmpbin[i] = pars[2]*ntmpbin[i]+pars[3]*ntmpbin2[i]+pars[4];
      }
      return ntmpbin;
   }

   SinglePulseShapeFunctor::SinglePulseShapeFunctor(const HcalPulseShapes::Shape& pulse) {
      for(int i=0;i<256;i++) {
         pulse_hist[i] = pulse(i);
      }
   }
  
   SinglePulseShapeFunctor::~SinglePulseShapeFunctor() {
   }
  
   double SinglePulseShapeFunctor::operator()(const std::vector<double>& pars) const {
      constexpr unsigned nbins = 10;
      unsigned i;

      //calculate chisquare
      double chisq = 0;
      double delta;
      std::array<float,nbins> pulse_shape = std::move(funcHPDShape(pars,pulse_hist));
      for (i=0;i<nbins; ++i) {
         delta = (psFit_y[i]- pulse_shape[i])/psFit_erry[i];
         chisq += delta*delta;
      }
      return chisq;
   }
  
   DoublePulseShapeFunctor::DoublePulseShapeFunctor(const HcalPulseShapes::Shape& pulse) {
      for(int i=0;i<256;i++) {
         pulse_hist[i] = pulse(i);
      }
   }

   DoublePulseShapeFunctor::~DoublePulseShapeFunctor() {
   }

   double DoublePulseShapeFunctor::operator()(const std::vector<double>& pars) const {
      constexpr unsigned nbins = 10;
      unsigned i;

      //calculate chisquare
      double chisq = 0;
      double delta;
      //double val[1];
      std::array<float,nbins> pulse_shape = std::move(func_DoublePulse_HPDShape(pars,pulse_hist));
      for (i=0;i<nbins; ++i) {
         delta = (psFit_y[i]- pulse_shape[i])/psFit_erry[i];
         chisq += delta*delta;
      }
      return chisq;
   }

   std::auto_ptr<SinglePulseShapeFunctor> spsfPtr_;
   std::auto_ptr<DoublePulseShapeFunctor> dpsfPtr_;

   double singlePulseShapeFunc( double x[3] ) {
      std::vector<double> pars;
      for(int i=0; i<3; i++) pars.push_back(x[i]);
      return (*spsfPtr_)(pars);
   }

   double doublePulseShapeFunc( double x[5] ) {
      std::vector<double> pars;
      for(int i=0; i<5; i++) pars.push_back(x[i]);
      return (*dpsfPtr_)(pars);
   }

}

PulseShapeFitOOTPileupCorrection::PulseShapeFitOOTPileupCorrection()
{
}

void PulseShapeFitOOTPileupCorrection::setPulseShapeTemplate(const HcalPulseShapes::Shape& ps) {
   spsf_.reset(new FitterFuncs::SinglePulseShapeFunctor(ps));
   dpsf_.reset(new FitterFuncs::DoublePulseShapeFunctor(ps));
   FitterFuncs::spsfPtr_.reset(new FitterFuncs::SinglePulseShapeFunctor(ps));
   FitterFuncs::dpsfPtr_.reset(new FitterFuncs::DoublePulseShapeFunctor(ps));
}

void PulseShapeFitOOTPileupCorrection::apply(const CaloSamples & cs, const std::vector<int> & capidvec,
                       const HcalCalibrations & calibs, std::vector<double> & correctedOutput) const
{
   FitterFuncs::cntNANinfit = 0;

//   CaloSamples cs;
//   coder.adc2fC(digi,cs);
   std::vector<double> chargeVec, pedVec;
   std::vector<double> energyVec, pedenVec;
   double TSTOT = 0, TStrig = 0; // in fC
   double TSTOTen = 0; // in GeV
   for(int ip=0; ip<cs.size(); ip++){
//      const int capid = digi[ip].capid();
      const int capid = capidvec[ip];
      double charge = cs[ip];
      double ped = calibs.pedestal(capid);
      double gain = calibs.respcorrgain(capid);

      double energy = charge*gain;
      double peden = ped*gain;

      chargeVec.push_back(charge); pedVec.push_back(ped);
      energyVec.push_back(energy); pedenVec.push_back(peden);

      TSTOT += charge - ped;
      TSTOTen += energy - peden;
      if( ip ==4 ){
         TStrig = charge - ped;
      }
   }
   std::vector<double> fitParsVec;
   if( TStrig >= 4 && TSTOT >= 10 ){
      pulseShapeFit(energyVec, pedenVec, chargeVec, pedVec, TSTOTen, fitParsVec, spsf_, dpsf_);
//      double time = fitParsVec[1], ampl = fitParsVec[0], uncorr_ampl = fitParsVec[0];
   }
   correctedOutput.swap(fitParsVec); correctedOutput.push_back(FitterFuncs::cntNANinfit);
}

int PulseShapeFitOOTPileupCorrection::pulseShapeFit(const std::vector<double> & energyVec, const std::vector<double> & pedenVec, const std::vector<double> &chargeVec, const std::vector<double> &pedVec, const double TSTOTen, std::vector<double> &fitParsVec, const std::auto_ptr<FitterFuncs::SinglePulseShapeFunctor>& spsf, const std::auto_ptr<FitterFuncs::DoublePulseShapeFunctor>& dpsf) const{

   int n_max=0;
   int n_above_thr=0;
   int first_above_thr_index=-1;
   int max_index[10]={0,0,0,0,0,0,0,0,0,0};

   double TSMAX=0;
   double TSMAX_NOPED=0;
   int i_tsmax=0;

   for(int i=0;i<10;i++){
      if(energyVec[i]>TSMAX){
         TSMAX=energyVec[i];
         TSMAX_NOPED=energyVec[i]-pedenVec[i];
         i_tsmax = i;
      }
   }

   double TIMES[10]={-100,-75,-50,-25,0,25,50,75,100,125};

   if(n_max==0){
      max_index[0]=i_tsmax;
   }

   double error = 1.;
   for(int i=0;i<10;i++){
      FitterFuncs::psFit_x[i]=i;
      FitterFuncs::psFit_y[i]=energyVec[i];
      FitterFuncs::psFit_erry[i]=error;
   }

   TFitterMinuit * gMinuit = new TFitterMinuit();
   gMinuit->SetPrintLevel(-1);

   for(int i=0;i!=10;++i){
      if((chargeVec[i])>6){
         n_above_thr++;
         if(first_above_thr_index==-1){
            first_above_thr_index=i;
         }
      }
   }

   // Fixed Maximum Finder
   for( int i=0 ; i < 10; ++i ) {
      switch( i ) {
         case 0:
            if(chargeVec[i]<=6 && chargeVec[i+1]<=6) continue;
            if( chargeVec[i+1] < chargeVec[i] ) {
               max_index[n_max++] = i;
            }
            break;
         case 9:
            if(chargeVec[i]<=6 && chargeVec[i-1]<=6) continue;
            if( chargeVec[i-1] < chargeVec[i] ) {
               max_index[n_max++] = i;
            }
            break;
         default:
            if(chargeVec[i-1]<=6 && chargeVec[i]<=6 && chargeVec[i+1]<=6) continue;
            if( chargeVec[i-1] < chargeVec[i] && chargeVec[i+1] < chargeVec[i]) {
               max_index[n_max++] = i;
            }
            break;
         }
      }

      if(n_max==0){
         max_index[0]=i_tsmax;
        //n_max=1; // there's still one max if you didn't find any...
      }

      double xmin_sp[3] = {}, xmin_dp[5] = {};
      double reqmin = 1.0E-12;

      int icount =0;
      int ifault =0;
      int kcount = 2500;
      int konvge = 10;
      int numres = 0;
      double ynewlo = 0;

      if(n_above_thr<=5){
         // Set starting values and step sizes for parameters
         double vstart[3] = {TIMES[i_tsmax-1],TSMAX_NOPED,0};

         double step[3] = {0.1,0.1,0.1};

         ynewlo = FitterFuncs::singlePulseShapeFunc(vstart);

         nelmin(FitterFuncs::singlePulseShapeFunc, 3, vstart, xmin_sp, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

      } else {
         if(n_max==1){
            // Set starting values and step sizes for parameters
            double vstart[5] = {TIMES[i_tsmax-1],TIMES[first_above_thr_index-1],TSMAX_NOPED,0,0};

            Double_t step[5] = {0.1,0.1,0.1,0.1,0.1};

            ynewlo = FitterFuncs::doublePulseShapeFunc(vstart);

            nelmin(FitterFuncs::doublePulseShapeFunc, 5, vstart, xmin_dp, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

        } else if(n_max>=2) {
           // Set starting values and step sizes for parameters
           double vstart[5] = {TIMES[max_index[0]-1],TIMES[max_index[1]-1],TSMAX_NOPED,0,0};

           double step[5] = {0.1,0.1,0.1,0.1,0.1};

           ynewlo = FitterFuncs::doublePulseShapeFunc(vstart);

           nelmin(FitterFuncs::doublePulseShapeFunc, 5, vstart, xmin_dp, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
        }
     }

     double timeval1fit=-999;
     double chargeval1fit=-999;
     double timeval2fit=-999;
     double chargeval2fit=-999;
     double pedvalfit=-999;

     if(n_above_thr<=5) {
        timeval1fit = xmin_sp[0];
        chargeval1fit = xmin_sp[1];
        pedvalfit = xmin_sp[2];
     } else {
        timeval1fit = xmin_dp[0];
        timeval2fit = xmin_dp[1];
        chargeval1fit = xmin_dp[2];
        chargeval2fit = xmin_dp[3];
        pedvalfit = xmin_dp[4];
     }

     int fitStatus = ifault;

     double timevalfit=0.;
     double chargevalfit=0.;
     if(n_above_thr<=5) {
        timevalfit=timeval1fit;
        chargevalfit=chargeval1fit;
     } else {
        if(fabs(timeval1fit)<fabs(timeval2fit)) {// if timeval1fit and timeval2fit are differnt, choose the one which is closer to zero
           timevalfit=timeval1fit;
           chargevalfit=chargeval1fit;
        } else if(fabs(timeval2fit)<fabs(timeval1fit)) {// if timeval1fit and timeval2fit are differnt, choose the one which is closer to zero
           timevalfit=timeval2fit;
           chargevalfit=chargeval2fit;
        } else if(timeval1fit==timeval2fit) { // if the two times are the same, then for charge we just sum the two  
           timevalfit=(timeval1fit+timeval2fit)/2;
           chargevalfit=chargeval1fit+chargeval2fit;
        } else {
           timevalfit=-999.;
           chargevalfit=-999.;
        }
     }

     fitParsVec.clear();

     fitParsVec.push_back(chargevalfit);
     fitParsVec.push_back(timevalfit);
     fitParsVec.push_back(pedvalfit);
     fitParsVec.push_back(ynewlo);
     fitParsVec.push_back(icount);
     fitParsVec.push_back(numres);
     fitParsVec.push_back(ifault);

     return fitStatus;
}
