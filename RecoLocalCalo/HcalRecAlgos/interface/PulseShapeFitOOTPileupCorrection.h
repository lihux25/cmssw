#ifndef PulseShapeFitOOTPileupCorrection_h
#define PulseShapeFitOOTPileupCorrection_h 1

#include <typeinfo>

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

#include <TMinuit.h>
#include "TFitterMinuit.h"

#include <TH1F.h>
#include "Minuit2/FCNBase.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "RecoLocalCalo/HcalRecAlgos/src/HybridMinimizer.h"

namespace FitterFuncs{
  
   class PulseShapeFunctor {
      public:
         PulseShapeFunctor(const HcalPulseShapes::Shape& pulse);
         ~PulseShapeFunctor();
         double EvalSinglePulse(const std::vector<double>& pars) const;
         double EvalDoublePulse(const std::vector<double>& pars) const;
      private:
         std::array<float,256> pulse_hist;
         std::vector<float> acc25nsVec, diff25nsItvlVec;
         std::vector<float> accVarLenIdxZEROVec, diffVarItvlIdxZEROVec;
         std::vector<float> accVarLenIdxMinusOneVec, diffVarItvlIdxMinusOneVec;
   };
   
}

class PulseShapeFitOOTPileupCorrection
{
public:
    PulseShapeFitOOTPileupCorrection();
    ~PulseShapeFitOOTPileupCorrection();

    // Main correction application method to be implemented by
    // derived classes. Arguments are as follows:
    //
    //
    // Some of the input arguments may be ignored by derived classes.
    //
    void apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalCalibrations & calibs, std::vector<double> & correctedOutput) const;

    void setPulseShapeTemplate(const HcalPulseShapes::Shape& ps);

private:

    bool useDataPulseShape_;

    int pulseShapeFit(const std::vector<double> & energyVec, const std::vector<double> & pedenVec, const std::vector<double> &chargeVec, const std::vector<double> &pedVec, const double TSTOTen, std::vector<double> &fitParsVec) const;

    PSFitter::HybridMinimizer * hybridfitter;

    int cntsetPulseShape;
};

#endif // PulseShapeFitOOTPileupCorrection_h
