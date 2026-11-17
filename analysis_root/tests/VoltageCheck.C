#include "../TPC.C"
#include <chrono>
// root -l -b -q 'analysis_root/tests/VoltageCheck.C()'
// Synthetic LM scan: no installed MS map or external voltage source needed.
void VoltageCheck() {
    using namespace PersonalTPC;
    auto folder=fs::temp_directory_path()/("tpc-voltage-check-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder);
    struct Cleanup {fs::path path;~Cleanup(){fs::remove_all(path);}} cleanup{folder};
    auto write=[](const fs::path& path,const std::string& text){std::ofstream f(path);f<<text;};
    auto db=folder/"cuts.ini",study=folder/"study.ini";
    write(db,"[global]\nsetup=large_medium\ndrift_voltage_v=-10.25\nreject=Time >= 228\n"
             "[run 1003]\ndrift_voltage_v=-10.75\n[run 1004]\ndrift_voltage_v=0\n"
             "[run 1005]\ndrift_voltage_v=0.123456789012345\n"
             "[run 1006]\nsetup=medium_small\ndrift_voltage_v=\n"
             "[run 1009]\ndrift_voltage_v=12.75\n"
             "[run 1605]\ndrift_voltage_v=\n[run 9999]\ndrift_voltage_v=999\n");
    write(study,"[global]\nreject=Channel == 42\n[run 1009]\nreject=Charge > 0\n");
    Config config(folder.c_str(),(folder/"results").c_str(),"",db.c_str(),study.c_str());
    auto makeRun=[&](int run,int events,bool burst=false) {
        TFile file((folder/("ntuple"+std::to_string(run)+".root")).c_str(),"RECREATE");
        TTree tree("ntuple","Synthetic voltage selection");
        int event,channel,hit=1,time,charge=100;
        tree.Branch("Event",&event);tree.Branch("Channel",&channel);tree.Branch("Hit",&hit);tree.Branch("Time",&time);tree.Branch("Charge",&charge);
        for(event=1;event<=events;++event) {
            time=10;
            if(burst && event==2)for(channel=0;channel<config.moduleChannels;++channel)tree.Fill();
            else {channel=152;tree.Fill();channel=42;tree.Fill();time=228;channel=152;tree.Fill();}
        }
        tree.Write();file.Close();
    };
    makeRun(1001,2);makeRun(1002,3);makeRun(1003,1);makeRun(1004,1);makeRun(1005,1);
    makeRun(1007,3,true);makeRun(1008,3);makeRun(1009,1);makeRun(1603,7);
    // Invalid trees for unknown runs prove skipping occurs before ROOT reads
    // or missing-MS-map checks. Unrelated names/directories are not inputs.
    write(folder/"ntuple1006.root","unknown input deliberately not a ROOT tree");
    write(folder/"ntuple1605.root","unknown input deliberately not a ROOT tree");
    write(folder/"notes.root","");write(folder/"ntuple10000.root.bak","");
    fs::create_directory(folder/"ntuple9999.root");
    std::ostringstream diagnostics;
    auto voltages=ReadVoltages(config,diagnostics);
    Check(voltages.size()==9 && voltages.at(1001)==-10.25 && voltages.at(1003)==-10.75 && voltages.at(1004)==0,"Discovery/inheritance/signed fractional voltages failed");
    Check(voltages.at(1005)==CutDatabase::Number("0.123456789012345"),"Discovery truncated DB precision");
    Check(!voltages.count(1006) && !voltages.count(1605) && !voltages.count(9999),"Unknown or unavailable run entered scan");
    Check(diagnostics.str().find("Run 1006:")!=std::string::npos && diagnostics.str().find("Run 1605:")!=std::string::npos,"Unknown runs lack skip diagnostic");
    Config override=config;override.setupOverride="medium_small";override.runInfo.clear();
    Check(ReadVoltages(override,diagnostics)==voltages,"Setup override changed DB voltages or required a missing map");
    Voltage(config);
    auto out=config.output/"voltage";
    std::ifstream csv(out/"run_summary.csv");std::string line;std::getline(csv,line);
    std::map<int,std::vector<std::string>> rows;std::set<int> selected;
    while(std::getline(csv,line)) {
        auto row=CSV(line);int run=std::stoi(row.at(0));rows[run]=row;
        Check(CutDatabase::Number(row.at(1))==voltages.at(run),"Export truncated/suppressed signed voltage");
        if(row.at(11)=="True")selected.insert(run);
    }
    Check(rows.size()==9 && selected==std::set<int>({1002,1003,1004,1005,1009}),"Largest-event selection, noise tie-break, lower-run tie-break or excluded1603 changed");
    Check(rows.at(1002).at(10)=="3" && rows.at(1009).at(10)=="0","Independent global/run study cuts did not apply");
    Check(rows.at(1007).at(3)=="2" && rows.at(1007).at(7)=="1","Noise tie-break fixture did not classify the burst");
    TFile file((out/"time_histograms.root").c_str(),"READ");
    for(int run:selected) {
        auto h=dynamic_cast<TH1D*>(file.Get(("time_run"+std::to_string(run)).c_str()));
        Check(h && h->GetEntries()==CutDatabase::Number(rows.at(run).at(10)),"Selected histogram missing or study cut count changed");
        auto info=dynamic_cast<TNamed*>(file.Get(("run"+std::to_string(run)+"_drift_voltage_v").c_str()));
        Check(info && CutDatabase::Number(info->GetTitle())==voltages.at(run),"ROOT configuration lost voltage precision");
    }
    Check(!file.Get("time_run1603") && !file.Get("time_run1008"),"Excluded/unselected run entered overlays");
    file.Close();
    auto signedSummary=Analyze(1003,voltages.at(1003),config,config.Mapping(1003));
    Check(signedSummary.voltage==-10.75 && Label(signedSummary).find("-10.75 V")==0,"Analyze/legend truncated voltage");
    auto emptyDB=folder/"unknown.ini";write(emptyDB,"[global]\ndrift_voltage_v=\n");
    Config unknown(folder.c_str(),(folder/"unknown-results").c_str(),"",emptyDB.c_str(),study.c_str());
    Check(ReadVoltages(unknown,diagnostics).empty(),"Empty voltage inherited a fabricated default");
    bool failed=false;
    try{Voltage(unknown);}catch(const std::exception& e){failed=std::string(e.what()).find("known drift_voltage_v")!=std::string::npos;}
    Check(failed,"Entirely unknown scan did not explain why it cannot run");
    std::cout<<"PASS: native DB-only voltage discovery, signed doubles/precision, unknown skips, setup override independence, study cuts, repeat-voltage winners/ties and excluded1603."<<std::endl;
}
