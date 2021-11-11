#ifndef MDSim_h
#define MDSim_h
#include "../declarations.h"
#include "SimConfig.h"
#include "SimStepData.h"
#include "AdditionalData.h"

class MDSim
{
    public:
    // Members
    SimConfig Cfg;
    AdditionalData AddedData;
    RandGen* rng;

    // Functions
    void ForcesFunc(SimConfig&, SimStepData&, AdditionalData&, int SampleInd);
    void PrintFunc(SimConfig&, SimStepData&, AdditionalData&);
    bool CheckFeedbackFunc(SimConfig&, SimStepData&, AdditionalData&);
    void FeedbackFunc(SimConfig&, SimStepData&, AdditionalData&);    
    void RunSim();
    
    MDSim() = delete;
    MDSim(SimConfig Cfg, RandGen* rng);

    private:
    void RunSim(SimConfig Cfg, AdditionalData AddedData);
};
#endif