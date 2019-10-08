// -*- C++ -*-
//
// Package:    HLTAnalysis/TriggerAnalyzerRAWMiniAOD
// Class:      TriggerAnalyzerRAWMiniAOD
// 
/**\class TriggerAnalyzerRAWMiniAOD TriggerAnalyzerRAWMiniAOD.cc HLTAnalysis/TriggerAnalyzerRAWMiniAOD/plugins/TriggerAnalyzerRAWMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Thomas
//         Created:  Fri, 24 Mar 2017 04:09:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TriggerAnalyzerRAWMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet&);
      ~TriggerAnalyzerRAWMiniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  bool PassOfflineMuonSelection(const pat::Muon *mu, reco::Vertex::Point PV);
  bool PassOfflineElectronSelection(const pat::Electron * ele, reco::Vertex::Point PV);
  bool RecoHLTMatching(const edm::Event&,double recoeta, double recophi, std::string filtername, double dRmatching = 0.3);
  double VarStudied( const edm::Event& iEvent, double recoeta, double recophi,edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > varToken_,  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> candToken_,   bool  dividebyE, bool dividebyEt, double dRmatching =0.3);

      // ----------member data ---------------------------


  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trgresultsHLT2Token_;

  edm::EDGetTokenT<std::vector<pat::Jet> > jet_token;
  edm::EDGetTokenT<std::vector<pat::Muon> > muon_token;
  edm::EDGetTokenT<std::vector<pat::Electron> > electron_token;
  edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> et_Filter_Token_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> showershape_Filter_Token_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> dphi_Filter_Token_;

  edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > showershape_Var_Token_;
  edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > hovere_Var_Token_;
  edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > trackiso_Var_Token_;

  edm::Service<TFileService> fs;
  
  TH1F* h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_den;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_num;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_numl1;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_den;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_num;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_den;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_num;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_den;
  TH1F* h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_num;
  TH1F* h_ele35wptight_lastfilter_den;
  TH1F* h_ele35wptight_lastfilter_num;
  TH1F* h_sietaieta_HLT;
  TH1F* h_hoe_HLT;
  TH1F* h_trackiso_HLT;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerAnalyzerRAWMiniAOD::TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet& iConfig)

{
  trigobjectsMINIAODToken_ = consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("slimmedPatTrigger"));
  trigobjectsRAWToken_=consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD::HLT2"));  

  trgresultsORIGToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );
  trgresultsHLT2Token_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT2") );


  showershape_Var_Token_  = consumes<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > ( edm::InputTag("hltEgammaClusterShape","sigmaIEtaIEta5x5","HLT2") );
  hovere_Var_Token_  = consumes<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > ( edm::InputTag("hltEgammaHoverE","","HLT2") );
  trackiso_Var_Token_  = consumes<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > ( edm::InputTag("hltEgammaEleGsfTrackIso","","HLT2")  );

  et_Filter_Token_ = consumes<trigger::TriggerFilterObjectWithRefs> ( edm::InputTag("hltEG35L1SingleEGOrEtFilter","","HLT2") ) ;
  showershape_Filter_Token_ = consumes<trigger::TriggerFilterObjectWithRefs> ( edm::InputTag("hltEle35noerWPTightClusterShapeFilter","","HLT2") );
  dphi_Filter_Token_ = consumes<trigger::TriggerFilterObjectWithRefs> ( edm::InputTag("hltEle35noerWPTightGsfDphiFilter","","HLT2") );
  

  jet_token = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets") );
  muon_token = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );
  electron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons") );
  PV_token = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  

  //now do what ever initialization is needed
  //   usesResource("TFileService");

  h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_den= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_den","",50,0,500);
  h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_num= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_num","",50,0,500);
  h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_numl1= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_numl1","",50,0,500);
  h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_den= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_den","",101,0,1.01);
  h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_num= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_num","",101,0,1.01);
  h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_den= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_den","",5,0,5);
  h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_num= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_num","",5,0,5);
  h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_den= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_den","",100,0,100);
  h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_num= fs->make<TH1F>("h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_num","",100,0,100);
  h_ele35wptight_lastfilter_den= fs->make<TH1F>("h_ele35wptight_lastfilter_den","",20,0,100);
  h_ele35wptight_lastfilter_num= fs->make<TH1F>("h_ele35wptight_lastfilter_num","",20,0,100);

  h_sietaieta_HLT= fs->make<TH1F>("h_sietaieta_HLT","",100,0,0.05);
  h_hoe_HLT= fs->make<TH1F>("h_hoe_HLT","",100,0,0.2);
  h_trackiso_HLT= fs->make<TH1F>("h_trackiso_HLT","",100,0,0.5);


}


TriggerAnalyzerRAWMiniAOD::~TriggerAnalyzerRAWMiniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnalyzerRAWMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace reco;
   using namespace std;



   // ****************Part 1. Accessing some trigger information ************* 
   bool passHLT_IsoMu24(false);
   bool passHLT_Mu3_PFJet200DeepCSV_1p59(false), passHLT_Mu3_L1SingleJet180(false), passHLT_PFJet200DeepCSV_1p59(false);   

   //Accessing trigger bits:
   //This works in both RAW, AOD or MINIAOD 
   //Here we access the decision provided by the HLT (i.e. original trigger step). 
   edm::Handle<edm::TriggerResults> trigResults;
   iEvent.getByToken(trgresultsORIGToken_, trigResults);
   if( !trigResults.failedToGet() ) {
     int N_Triggers = trigResults->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResults.product()->accept(i_Trig)) {
	 TString TrigPath =trigName.triggerName(i_Trig);
	 //	 cout << "Passed path: " << TrigPath<<endl;
	 if(TrigPath.Index("HLT_IsoMu24_v") >=0) passHLT_IsoMu24=true; 
	 //Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu24_v"
       }
     }
   }
   //Exercise 1: 
   //Clone and *then* modify the code above in order to save the decision of your customized HLT menu in the booleans passHLT_Mu3_PFJet200DeepCSV_1p59, passHLT_Mu3_L1SingleJet180, passHLT_PFJet200DeepCSV_1p59
   //Do not directly edit the code above as you will also need the use the original HLT_IsoMu24 decision later on.

   



   //Accessing the trigger objects in MINIAOD
   //This recipe works for MINIAOD only
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

   const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
   for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
     obj.unpackFilterLabels(iEvent,*trigResults);
     obj.unpackPathNames(names);
     for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
       string myfillabl=obj.filterLabels()[h];
       //cout << "Trigger object name, pt, eta, phi: "
       //	    << myfillabl<<", " << obj.pt()<<", "<<obj.eta()<<", "<<obj.phi() << endl;
     }
   }

   //Exercise 2: uncomment the lines above to print all the trigger objects and their corresponding pt, eta, phi. 
   

   //Accessing the trigger objects in RAW/AOD
   //Printing here all trigger objects corresponding to the filter hltL3MuFiltered3
   edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
   iEvent.getByToken(trigobjectsRAWToken_ ,triggerObjectsSummary);
   trigger::TriggerObjectCollection selectedObjects;
   if (triggerObjectsSummary.isValid()) {
     size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltL3MuFiltered3","","HLT2") );
     trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
     if (filterIndex < (*triggerObjectsSummary).sizeFilters()) { 
       const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
       for (size_t j = 0; j < keys.size(); j++) {
	 //trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	 //cout <<"object found, printing pt, eta, phi: " <<foundObject.pt()<<", "<<foundObject.eta()<<", "<< foundObject.phi() <<endl;
       }
     }
   }
   
   //Exercise 3: uncomment the two lines above and modify the input tag to print all trigger objects corresponding to the last filter of the HLT_Mu3_PFJet200DeepCSV_1p59 path (btagged jet with pt>200 GeV)



   // **************** Part 2. Accessing some offline information ************** 
   
   //What you really want to do is to assess the trigger performances on top of an offline selection. 
    
   
   //Offline jets
   //Find the highest pt b jet (medium WP i.e. csv>0.8484 for b tagging) in the event.
   //Find the highest csv of a jet with pt>250 GeV in the event.
   //Count the nb of bjets with pt>200 GeV in the event
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByToken(jet_token,jets );

   double leadingbjetpt(-100), leadingbjeteta(-100),leadingbjetphi(-100); 
   double highestcsv_jetpt250 =-1;
   int nbjetspt200 = 0;
   for( std::vector<pat::Jet>::const_iterator jet = (*jets).begin(); jet != (*jets).end(); jet++ ) {
     double ptjet = jet->pt();
     double etajet = jet->eta();
     double phijet = jet->phi();
     double csvjet = jet->bDiscriminator("pfDeepCSVJetTags:probb")+ jet->bDiscriminator("pfDeepCSVJetTags:probbb");//cf https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
     //The next following lines just remove e/mu from (semi)leptonic ttbar. 
     if( jet->muonEnergyFraction() >0.7)continue;
     if( jet->electronEnergyFraction() >0.7)continue;
     
     if(abs( etajet)>2.4) continue; //Only consider jets in tracker acceptance since we want to do b tagging. 
     if(csvjet>0.4184&& ptjet>leadingbjetpt) { leadingbjetpt = ptjet; leadingbjeteta =etajet; leadingbjetphi = phijet;} 
     if(ptjet>250&& csvjet>highestcsv_jetpt250) { highestcsv_jetpt250 = csvjet;} 
     if(csvjet>0.4184&& ptjet>200) { nbjetspt200++; } 
   }

   
   //Offline muons 
   edm::Handle< std::vector<pat::Muon> > muons;
   iEvent.getByToken(muon_token,muons );
   //We also need the vertices here
   edm::Handle<std::vector<Vertex> > theVertices;
   iEvent.getByToken(PV_token,theVertices) ;
   int nvertex = theVertices->size();
   Vertex::Point PV(0,0,0);
   if( nvertex) PV = theVertices->begin()->position();
   //Count the nb of offline muons with pt >3
   //Find the highest pt muon
   int nmuonspt3 =0; 
   double leadingmuonpt(-10),leadingmuoneta(-10),leadingmuonphi(-10);
   for( std::vector<pat::Muon>::const_iterator muon = (*muons).begin(); muon != (*muons).end(); muon++ ) {
     if(!PassOfflineMuonSelection(&*muon,PV)) continue;
     double ptmuon = muon->pt();
     double etamuon = muon->eta();
     double phimuon = muon->phi();
     if(ptmuon>=3) nmuonspt3++;
     if(ptmuon>leadingmuonpt){leadingmuonpt=ptmuon;leadingmuoneta=etamuon;leadingmuonphi=phimuon;}
   }


   
   //We can now fill some histograms (numerator and denominator to study the efficiency of our favourite path).
   //Here we factorize the muon and jet legs and measure their efficiencies separately
    

   //Effcy vs pt of the leading b jet: 
   h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_den->Fill(leadingbjetpt);
   if(passHLT_Mu3_PFJet200DeepCSV_1p59) h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_num->Fill(leadingbjetpt);
   if(passHLT_Mu3_L1SingleJet180) h_mu3pfjet200deepcsv1p59_vs_leadbjetpt_numl1->Fill(leadingbjetpt); //For L1 turn on only

   //Effcy vs the highest csv of jet with pt>250:
   h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_den->Fill(highestcsv_jetpt250);
   if(passHLT_Mu3_PFJet200DeepCSV_1p59) h_mu3pfjet200deepcsv1p59_vs_highestcsv_jetpt250_num->Fill(highestcsv_jetpt250);


   //Effcy vs nb of bjets with pt>200 
   h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_den->Fill(nbjetspt200);
   if(passHLT_Mu3_PFJet200DeepCSV_1p59) h_mu3pfjet200deepcsv1p59_vs_nbjetspt200_num->Fill(nbjetspt200);

   //Effcy vs leading muon pt:
   h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_den->Fill(leadingmuonpt);
   if(passHLT_Mu3_PFJet200DeepCSV_1p59) h_mu3pfjet200deepcsv1p59_vs_leadingmuonpt_num->Fill(leadingmuonpt);

   //Exercise 4: Define the denominator (offline selection) for the histograms filled above and look at the obtained efficiency plots. 
   //For this exercise, you need to run on at least 10k events in order to start to see the jet turn on. 
   

   
   //Exercise 5: Take a look at the histograms obtained in the following root file, with much higher stat: 
   ///afs/cern.ch/work/l/lathomas/public/HLTTutorial_27Oct2017/OutputFiles/out_singlemuon_highstat.root
   //Check the efficiency vs the leading jet pt. 

   //Questions: 
   //- Would you say that the plateau efficiency represents the probability for a b quark jet to fire the b tagging condition at HLT? 
   //- Do you understand why the efficiency increases with the n(bjets)? 
   //-Is the efficiency measured vs muon pt unbiased when running on a SingleMuon dataset?  


   
   

   //Let's finally see a Tag and Probe example on Z(ee)
   //Offline electrons
   edm::Handle< std::vector<pat::Electron> > electrons;
   iEvent.getByToken(electron_token,electrons );
   //First loop to find a tag electron
   for( std::vector<pat::Electron>::const_iterator tagele = (*electrons).begin(); tagele != (*electrons).end(); tagele++ ) {
     if(!PassOfflineElectronSelection(&*tagele,PV)) continue;
     double pttagele = tagele->pt();
     double etatagele = tagele->eta();
     double phitagele = tagele->phi();
     if(pttagele<35) continue; 
     if(!RecoHLTMatching(iEvent,etatagele,phitagele,"hltEle35noerWPTightGsfTrackIsoFilter") ) continue;
     //We want to match the tag to the last filter of the HLT_Ele35_WPTight_Gsf
     //Take a look at  HLTrigger/Configuration/python/HLT_TutoEle35WPTight_cff.py to confirm that "hltEle35noerWPTightGsfTrackIsoFilter" is indeed the last filter of that path
     
     //Second loop on the probe
     for( std::vector<pat::Electron>::const_iterator probeele = (*electrons).begin(); probeele != (*electrons).end(); probeele++ ) {
       if(tagele==probeele)continue;//Tag and Probe should be different (obviously)
       if(!PassOfflineElectronSelection(&*probeele,PV)) continue;
       double ptprobeele = probeele->pt();
       double etaprobeele = probeele->eta();
       double phiprobeele = probeele->phi();

       TLorentzVector p4tag, p4probe; 
       p4tag.SetPtEtaPhiM(pttagele,etatagele,phitagele,0);
       p4probe.SetPtEtaPhiM(ptprobeele,etaprobeele,phiprobeele,0);
       double mass = (p4tag+p4probe).Mag();
       
       if(mass<60 ||mass>120) continue;
       if(tagele->charge() * probeele->charge()>0)continue;
       h_ele35wptight_lastfilter_den->Fill(ptprobeele);
       if(RecoHLTMatching(iEvent,etaprobeele,phiprobeele,"hltEle35noerWPTightGsfTrackIsoFilter") ) h_ele35wptight_lastfilter_num->Fill(ptprobeele);

       
       //Note that everything above is done on miniAOD.
       //If you want to do something similar for a new path, then you need to rerun HLT from RAW (obviously) and modify a bit the RecoHLTMatching function when you retrieve the trigger objects. 
       //The previous examples should make it clear on how to do that. 

     
       
       //Now, new (advanced) topic: access the HLT ID variables for an electron/photon.
       //These are typically not stored, so you have to rerun HLT. 
       //The various filters in the Single electron triggers are (in that order): Et cut, sigmaietaieta cut, H/E cut, ..., dphi cut, track iso cut
       //We will take advantage that all these cuts are applied using the same EDFilter http://cmslxr.fnal.gov/source/HLTrigger/Egamma/src/HLTGeneric(QuadraticEta)Filter.cc 
       //Here we will check the distribution of each variable just before the cut is applied. 
       //For that we will need to specify the variable name and the previous filter name. All the work is done in this VarStudied function which essentially copies what is done in HLTGeneric(QuadraticEta)Filter
       
       
       bool  dividebyE = false; bool dividebyEt = false; 
       
       //First sietaieta
       double sietaieta_HLT = VarStudied(iEvent, etaprobeele,phiprobeele,showershape_Var_Token_,et_Filter_Token_ ,dividebyE, dividebyEt);
       h_sietaieta_HLT->Fill(sietaieta_HLT);
       //Next: H/E:
       dividebyE = true; dividebyEt = false; //For H/E
       double hoe_HLT  = VarStudied(iEvent, etaprobeele,phiprobeele,hovere_Var_Token_,showershape_Filter_Token_ ,dividebyE, dividebyEt);
       h_hoe_HLT->Fill(hoe_HLT);
       //Finally: trackiso 
       dividebyE =false; dividebyEt = true; //For isolation 
       double trackiso_HLT  = VarStudied(iEvent, etaprobeele,phiprobeele,trackiso_Var_Token_,dphi_Filter_Token_ ,dividebyE, dividebyEt);
       h_trackiso_HLT->Fill(trackiso_HLT);

     }
   }

   //Exercise 6.
   //Rerun the HLT_Ele35_WPTight_Gsf and to run on the SingleElectron dataset. 
   //In HLTrigger/Configuration/test do: 
   //cp HLT2_HLT.py HLT2_HLT_SingleEle.py
   //Open HLT2_HLT_SingleEle.py and make the following changes: 
   //The line:
   // process.load('HLTrigger.Configuration.HLT_TutoEffcySession_cff') 
   // should be replaced by: 
   // process.load('HLTrigger.Configuration.HLT_TutoEle35WPTight_cff')
   //The input files should be updated: 
   //fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/MINIAOD/22Jan2019-v2/70001/D8790A22-9BE2-624F-A1CC-6A4A7CAC82D7.root')
   //secondaryFileNames = cms.untracked.vstring(
   //                                         'root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/323/755/00000/F9080D12-2CDA-4F43-B7A9-C2F1D05E9C10.root', 
   //                                         'root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/323/755/00000/72D435E1-78D4-8047-92BD-3553117A5594.root',                       
   //                                         'root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/323/755/00000/266F4E2C-B30A-724E-9183-AD9368A10D51.root' 
   //)     

   //Make the above distributions (e.g. for sietaieta) only for probes passing the triggers and observe that the distributions are indeed truncated at the cut values reported in the HLT config file (HLTrigger/Configuration/python/HLT_TutoEle35WPTight_cff.py

   


}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnalyzerRAWMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool TriggerAnalyzerRAWMiniAOD::PassOfflineMuonSelection(const pat::Muon *mu, reco::Vertex::Point PV){
  if ( !(mu->isGlobalMuon() || mu->isTrackerMuon() )) return false;
  if ( !(mu->isPFMuon()) ) return false;
  const reco::TrackRef innerTrack = mu->innerTrack();
  if( innerTrack.isNull() )return false;
  
  bool goodGlb =  mu->isGlobalMuon() &&  mu->globalTrack()->normalizedChi2() < 3
    &&  mu->combinedQuality().chi2LocalPosition < 12  && mu->combinedQuality().trkKink < 20;
  bool good =  mu->innerTrack()->validFraction() >= 0.8 &&  mu->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451)  ;
  
  if(!good) return false;
  if(TMath::Abs(innerTrack->dxy(PV)) >0.1 ) return false;
  if(TMath::Abs(innerTrack->dz(PV)) >0.1 ) return false;  
  
  double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
  double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
  double photonIso = mu->pfIsolationR03().sumPhotonEt;
  
  double beta = mu->pfIsolationR03().sumPUPt;
  double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
  
  if(pfRelIsoMu >0.4) return false;
  return true;
}









bool TriggerAnalyzerRAWMiniAOD::PassOfflineElectronSelection(const pat::Electron * ele, reco::Vertex::Point PV){
  const reco::GsfTrackRef gsfTrack = ele->gsfTrack();
  if (!gsfTrack.isNonnull()) return false;
  if( TMath::Abs(gsfTrack->dxy(PV)) > 0.05  )  return false;
  if( TMath::Abs(gsfTrack->dz(PV)) > 0.2  )  return false;
  if(TMath::Abs(ele->superCluster()->eta()) >2.5) return false;
  else if( TMath::Abs(ele->superCluster()->eta()) < 1.479  ) {
    if( TMath::Abs(ele->full5x5_sigmaIetaIeta()) > 0.0103 ) return  false;
    if( TMath::Abs(ele->deltaEtaSuperClusterTrackAtVtx()) > 0.0105 ) return  false;
    if( TMath::Abs(ele->deltaPhiSuperClusterTrackAtVtx()) > 0.115  ) return  false;
    if( TMath::Abs(ele->hadronicOverEm())  > 0.104  ) return  false;
    if( TMath::Abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy() )>0.102 ) return  false;
    if( TMath::Abs(gsfTrack->dxy(PV)) > 0.0261) return false;
    if( TMath::Abs(gsfTrack->dz(PV)) > 0.41) return false;
    if(gsfTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>2) return  false;
  }
  else {
    if( TMath::Abs(ele->full5x5_sigmaIetaIeta()) > 0.0301 ) return  false;
    if( TMath::Abs(ele->deltaEtaSuperClusterTrackAtVtx()) > 0.00814 ) return  false;
    if( TMath::Abs(ele->deltaPhiSuperClusterTrackAtVtx()) > 0.182  ) return  false;
    if( TMath::Abs(ele->hadronicOverEm())  > 0.0897  ) return  false;
    if( TMath::Abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy() )>0.126 ) return  false;
    if( TMath::Abs(gsfTrack->dxy(PV)) > 0.0118) return false;
    if( TMath::Abs(gsfTrack->dz(PV)) > 0.822) return false;
    if(gsfTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>1) return  false;


  }

  double iso = (ele->pfIsolationVariables().sumChargedHadronPt 
		+ TMath::Max(0.0, ele->pfIsolationVariables().sumNeutralHadronEt + ele->pfIsolationVariables().sumPhotonEt - 0.5*ele->pfIsolationVariables().sumPUPt ) 
		) /ele->pt() ; 
  if(iso>0.2) return false; 
  return true;



}


bool TriggerAnalyzerRAWMiniAOD::RecoHLTMatching(const edm::Event& iEvent, double recoeta, double recophi, std::string filtername, double dRmatching){
  //In the next few lines one loops over all the trigger objects (corresponding to a given filter) and check whether one of them matches the reco object under study                                       
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

  const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackFilterLabels(iEvent,*trigResults);
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      std::string myfillabl=obj.filterLabels()[h];
      if( myfillabl.find(filtername)!=std::string::npos   && deltaR(recoeta,recophi, obj.eta(),obj.phi())<dRmatching ) return true;
    }
  }

  return false;
}




double TriggerAnalyzerRAWMiniAOD::VarStudied( const edm::Event& iEvent, double recoeta, double recophi,
					      edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > varToken_,  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> candToken_,   bool  dividebyE, bool dividebyEt, double dRmatching ){

  double thevar = 0.;

  //Inspired from http://cmslxr.fnal.gov/source/HLTrigger/Egamma/src/HLTGenericFilter.cc        
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  iEvent.getByToken (candToken_, PrevFilterOutput);

  edm::Handle<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > depMap;
  iEvent.getByToken (varToken_,depMap);

  std::vector<edm::Ref<std::vector<reco::RecoEcalCandidate> > > recoCands;


  if(PrevFilterOutput.isValid()&&  depMap.isValid() ){

    PrevFilterOutput->getObjects(trigger::TriggerCluster, recoCands);
    if(recoCands.empty())PrevFilterOutput->getObjects(trigger::TriggerPhoton, recoCands);

    double dRmin = dRmatching;
    for (unsigned int i=0; i<recoCands.size(); i++) {
      edm::Ref<std::vector<reco::RecoEcalCandidate> > ref = recoCands[i];
      typename edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > >::const_iterator mapi = (*depMap).find( ref );
      float vali = mapi->val;
      float EtaSC = ref->eta();
      float PhiSC = ref->phi();

      if(deltaR(recoeta,recophi,EtaSC,PhiSC ) > dRmin )continue;
      dRmin = deltaR(recoeta,recophi,EtaSC,PhiSC ) ;
      float energy = ref->superCluster()->energy();
      float et = ref->superCluster()->energy() * sin (2*atan(exp(-ref->eta())));
      thevar = (double) vali;
      if(dividebyE)thevar = (double)vali/energy;
      if(dividebyEt)thevar =(double) vali/et;

    }
  }
  return thevar;
}



//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
