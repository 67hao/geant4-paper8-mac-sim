#include "MyRunAction.hh"
#include "MyDetectorConstruction.hh"   // to read current sample material name
#include "MyPrimaryGeneratorAction.hh" // used directly below (static_cast<const MyPrimaryGeneratorAction*>)
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

// to fetch DetPlane's outer radius -> surface area for the crossing-flux tally
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"

#include <cmath>
#include <ctime>
#include <fstream>
#include <atomic>
#include <mutex>
#include <cstdio>
#include <iomanip>

namespace {
    std::atomic<G4double> g_primaryEnergy{ 0.0 };

    std::mutex g_materialMutex;
    G4String g_sampleMaterial = "S1";

    void SetGlobalSampleMaterial(const G4String& name) {
        std::lock_guard<std::mutex> lock(g_materialMutex);
        g_sampleMaterial = name;
    }
    G4String GetGlobalSampleMaterial() {
        std::lock_guard<std::mutex> lock(g_materialMutex);
        return g_sampleMaterial;
    }

    std::atomic<G4double> g_sampleThickness{ 0.0 };
    void SetGlobalSampleThickness(G4double t) {
        g_sampleThickness.store(t, std::memory_order_relaxed);
    }
    G4double GetGlobalSampleThickness() {
        return g_sampleThickness.load(std::memory_order_relaxed);
    }
}

// =========================================================================
// Constructor - NEW: khong con dang ky accumulable Edep/Pulse, khong con
// CreateH1/CreateNtuple cho Edep/Pulse -> giam overhead moi hit.
// NEW: khong con dang ky accumulable track-length (fTrackLengthNaI /
// fTrackLengthUncollidedNaI) - logic track-length tally da bi BO HOAN TOAN,
// chi con lai surface-crossing flux tally (T-Cross-like) tai DetPlane.
// =========================================================================
MyRunAction::MyRunAction()
    : G4UserRunAction(),
    fNCompton(0.0), fNPhotoelectric(0.0), fNPairProduction(0.0),
    fNRayleigh(0.0), fNOther(0.0), fNTotal(0.0)
{
    auto* accMan = G4AccumulableManager::Instance();
    accMan->RegisterAccumulable(fNCompton);
    accMan->RegisterAccumulable(fNPhotoelectric);
    accMan->RegisterAccumulable(fNPairProduction);
    accMan->RegisterAccumulable(fNRayleigh);
    accMan->RegisterAccumulable(fNOther);
    accMan->RegisterAccumulable(fNTotal);

    // surface-crossing (T-Cross-like) tally at DetPlane
    accMan->RegisterAccumulable(fFluxCrossNaI);
    accMan->RegisterAccumulable(fFluxCrossUncollidedNaI);

    // AnalysisManager: ten file .root dat theo macro /analysis/setFileName
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);
    analysisManager->SetDefaultFileType("root");

    // --- ROOT output ---
    // H1 0: pulse-height spectrum in NaI (1 keV bins, covers E0 up to 15 MeV)
    // H1 1: photon spectrum crossing DetPlane, weighted by flux = weight/|cosTheta|
    analysisManager->CreateH1("Edep_NaI", "Deposited energy in NaI", 15500, 0., 15500., "keV");
    analysisManager->CreateH1("E_cross", "Energy at DetPlane crossing (flux-weighted)", 15500, 0., 15500., "keV");

    // Ntuple 0 "Hits": same name/column as Paper2/3, so Count.cc works unchanged
    analysisManager->CreateNtuple("Hits", "Energy deposit in NaI per event");
    analysisManager->CreateNtupleDColumn("Edep_keV");
    analysisManager->FinishNtuple();

    // Ntuple 1 "Cross": one row per particle crossing DetPlane
    analysisManager->CreateNtuple("Cross", "DetPlane crossings (T-Cross-like)");
    analysisManager->CreateNtupleDColumn("E_keV");
    analysisManager->CreateNtupleDColumn("CosTheta");
    analysisManager->CreateNtupleDColumn("Flux");          // weight/|cosTheta|
    analysisManager->CreateNtupleIColumn("Uncollided");    // 1 if |E - E0| <= window*E0
    analysisManager->FinishNtuple();
}

MyRunAction::~MyRunAction() = default;

void MyRunAction::BeginOfRunAction(const G4Run*) {
    G4AccumulableManager::Instance()->Reset();
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Reset();

    auto detector = static_cast<const MyDetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    G4String matName = detector ? detector->GetSampleMaterialName() : "Unknown";
    SetGlobalSampleMaterial(matName);

    G4double thicknessMM = detector ? detector->GetSampleThickness() / mm : 0.0;
    SetGlobalSampleThickness(thicknessMM);

    auto pga = static_cast<const MyPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    if (pga) {
        g_primaryEnergy.store(pga->GetParticleEnergy(), std::memory_order_relaxed);
    }
    G4double energyKeV = g_primaryEnergy.load(std::memory_order_relaxed) / keV;

    G4String presetName = analysisManager->GetFileName();
    if (!presetName.empty()) {
        analysisManager->OpenFile();
    }
    else {
        char tag[128];
        std::snprintf(tag, sizeof(tag), "_t%.2fmm_E%.1fkeV_output", thicknessMM, energyKeV);
        analysisManager->OpenFile(matName + tag);
    }

    if (IsMaster()) {
        fTimer.Start();
        std::time_t t = std::time(nullptr);
        char buf[64];
        std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
        G4cout << "\n========== Run started: " << buf << " ==========" << G4endl;
    }
}

// =========================================================================
// EndOfRunAction - NEW: khong con edep/pulse, khong con track-length tally
// (F4-like) trong NaI. Chi con lai surface-crossing flux tally (T-Cross-like)
// tai DetPlane + process counting o Sample.
// =========================================================================
void MyRunAction::EndOfRunAction(const G4Run* run) {
    G4AccumulableManager::Instance()->Merge();
    auto* analysisManager = G4AnalysisManager::Instance();

    if (run->GetNumberOfEvent() == 0) return;

    G4int    nEvents = run->GetNumberOfEvent();
    G4double nCompt = fNCompton.GetValue();
    G4double nPhot = fNPhotoelectric.GetValue();
    G4double nPair = fNPairProduction.GetValue();
    G4double nRayl = fNRayleigh.GetValue();
    G4double nOther = fNOther.GetValue();
    G4double nTot = fNTotal.GetValue();

    // =====================================================================
    // Surface-crossing flux tally (T-Cross-like) at DetPlane, dat sat truoc
    // mat NaI. Chia cho DIEN TICH mat cat cua DetPlane (dung ban kinh ngoai
    // cua G4Tubs "DetPlane" duoc dat trong MyDetectorConstruction::Construct()).
    // =====================================================================
    G4double crossFluxRaw            = fFluxCrossNaI.GetValue();
    G4double crossFluxUncollidedRaw  = fFluxCrossUncollidedNaI.GetValue();

    G4double detPlaneArea_cm2 = 0.0;
    auto detPlaneLV = G4LogicalVolumeStore::GetInstance()->GetVolume("DetPlane_LV", false);
    if (detPlaneLV && detPlaneLV->GetSolid()) {
        if (auto tubs = dynamic_cast<G4Tubs*>(detPlaneLV->GetSolid())) {
            G4double rOuter = tubs->GetOuterRadius();
            detPlaneArea_cm2 = (CLHEP::pi * rOuter * rOuter) / (cm * cm);
        }
    }

    G4double crossFluxTotal      = 0.0; // /cm^2/particle - tuong duong output "flux" cua PHITS [T-Cross]
    G4double crossFluxUncollided = 0.0;
    if (detPlaneArea_cm2 > 0. && nEvents > 0) {
        crossFluxTotal      = crossFluxRaw           / detPlaneArea_cm2 / nEvents;
        crossFluxUncollided = crossFluxUncollidedRaw  / detPlaneArea_cm2 / nEvents;
    }

    if (IsMaster()) {
        G4cout << "\n================ ANALYSIS RESULTS ================" << G4endl;
        G4cout << "Sample material  : " << GetGlobalSampleMaterial() << G4endl;
        G4cout << "Number of events: " << nEvents << G4endl;

        G4cout << "\n--- Surface-Crossing Flux Tally (T-Cross-like) at DetPlane ---" << G4endl;
        G4cout << "  DetPlane area          = " << detPlaneArea_cm2 << " cm^2" << G4endl;
        G4cout << "  Flux total (crossing)  = " << crossFluxTotal << " /cm^2/particle" << G4endl;
        G4cout << "  Flux uncollided        = " << crossFluxUncollided
               << " /cm^2/particle  (window = +/-" << fEnergyWindowFrac * 100.0 << "% around E0)" << G4endl;

        G4cout << "\n--- First Interaction in Sample ---" << G4endl;
        G4cout << "  Compton        = " << nCompt << G4endl;
        G4cout << "  Photoelectric  = " << nPhot << G4endl;
        G4cout << "  Pair Production = " << nPair << G4endl;
        G4cout << "  Rayleigh       = " << nRayl << G4endl;
        G4cout << "  Other          = " << nOther << G4endl;
        G4cout << "  Total          = " << nTot << G4endl;

        if (nTot > 0) {
            G4cout << "\n--- R = N_process / (N_total - N_Rayleigh) ---" << G4endl;
            G4cout << "  R_Compton        = " << nCompt / (nTot - nRayl) << G4endl;
            G4cout << "  R_Photoelectric  = " << nPhot / (nTot - nRayl) << G4endl;
            G4cout << "  R_PairProduction = " << nPair / (nTot - nRayl) << G4endl;
            G4cout << "  R_Rayleigh       = " << nRayl / (nTot - nRayl) << G4endl;
            G4cout << "  Ratio Checking   = " << (nCompt + nPhot + nPair + nRayl + nOther) / nTot << G4endl;
        }
        G4cout << "==================================================\n" << G4endl;

        fTimer.Stop();
        std::time_t t = std::time(nullptr);
        char buf[64];
        std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
        G4double real = fTimer.GetRealElapsed();
        G4int h = (G4int)(real / 3600);
        G4int m = (G4int)((real - h * 3600) / 60);
        G4int s = (G4int)(real - h * 3600 - m * 60);
        G4cout << "========== Run ended  : " << buf << " ==========" << G4endl;
        G4cout << "========== Elapsed    : "
            << h << "h " << m << "m " << s << "s"
            << "  (real: " << std::fixed << std::setprecision(1) << real << " s)"
            << " ==========" << G4endl;

        G4double energy_keV = g_primaryEnergy.load(std::memory_order_relaxed) / keV;
        G4String matName = GetGlobalSampleMaterial();

        // --- CSV output: KHONG con Pulse_NaI, Edep_MeV, TrackLength_* ---
        std::string csvFile = "process_results.csv";

        std::ifstream checkFile(csvFile);
        bool writeHeader = !checkFile.good();
        checkFile.close();

        std::ofstream csv(csvFile, std::ios::app);
        if (writeHeader) {
            csv << "Material,Thickness_mm,Energy_keV,N_events,N_compton,N_photo,N_pair,N_rayleigh,N_other,N_total,"
                << "R_compton,R_photo,R_pair,R_rayleigh,"
                << "DetPlaneArea_cm2,CrossFlux_total,CrossFlux_uncollided" << std::endl;
        }

        G4double rC = (nTot > 0) ? nCompt / (nTot - nRayl) : 0;
        G4double rP = (nTot > 0) ? nPhot / (nTot - nRayl) : 0;
        G4double rPP = (nTot > 0) ? nPair / (nTot - nRayl) : 0;
        G4double rR = (nTot > 0) ? nRayl / (nTot - nRayl) : 0;

        G4double thicknessMM = GetGlobalSampleThickness();

        csv << matName << ","
            << std::fixed
            << std::setprecision(3) << thicknessMM << ","
            << std::setprecision(2)
            << energy_keV << ","
            << nEvents << ","
            << std::setprecision(0)
            << nCompt << "," << nPhot << "," << nPair << "," << nRayl << "," << nOther << "," << nTot << ","
            << std::setprecision(8)
            << rC << "," << rP << "," << rPP << "," << rR << ","
            << std::setprecision(6) << detPlaneArea_cm2 << ","
            << std::setprecision(8) << crossFluxTotal << ","
            << std::setprecision(8) << crossFluxUncollided
            << std::endl;
        csv.close();
    }

    analysisManager->Write();
    analysisManager->CloseFile();
}

void MyRunAction::AddProcessCount(const G4String& processName) {
    fNTotal += 1.0;
    if (processName == "compt") fNCompton += 1.0;
    else if (processName == "phot")  fNPhotoelectric += 1.0;
    else if (processName == "conv")  fNPairProduction += 1.0;
    else if (processName == "Rayl")  fNRayleigh += 1.0;
    else                             fNOther += 1.0;
}

// NEW: surface-crossing (T-Cross-like) tally accumulators
void MyRunAction::AddCrossingFluxNaIDet(G4double flux) {
    fFluxCrossNaI += flux;
}

void MyRunAction::AddCrossingFluxUncollidedNaIDet(G4double flux) {
    fFluxCrossUncollidedNaI += flux;
}

// NEW: ROOT output fillers (called from MyEventAction::EndOfEventAction)
void MyRunAction::FillEdepNaI(G4double edep) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(0, edep);
    analysisManager->FillNtupleDColumn(0, 0, edep / keV);
    analysisManager->AddNtupleRow(0);
}

void MyRunAction::FillCrossing(G4double ekin, G4double cosTheta,
                               G4double flux, G4bool uncollided) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(1, ekin, flux);
    analysisManager->FillNtupleDColumn(1, 0, ekin / keV);
    analysisManager->FillNtupleDColumn(1, 1, cosTheta);
    analysisManager->FillNtupleDColumn(1, 2, flux);
    analysisManager->FillNtupleIColumn(1, 3, uncollided ? 1 : 0);
    analysisManager->AddNtupleRow(1);
}

G4double MyRunAction::GetPrimaryEnergy() const {
    return g_primaryEnergy.load(std::memory_order_relaxed);
}
